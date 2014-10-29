function krk_output = readKnossosNml( krk_fname, krk_keepNodeAsStruct )
    
% READKNOSSOSNML: Read .nml files from Knossos into Matlab
%   
%   The function has the following arguments:
%       KRK_FNAME: Give the complete path of the file you want to read as a
%           string, e.g. 'E:\e_k0563\ribbons\skeleton_seeds\All_ek0563_Skeletons\ek0563-seed0001-khaase.002.nml'
%       KRK_KEEPNODEASSTRUCT: Optional! Standard version is 1. The numbers 0 or 1 indicate
%           whether the nodes should be saved as a struct additionally.
%       
%   => readKnossosNml( 'E:\e_k0563\ribbons\skeleton_seeds\All_ek0563_Skeletons\ek0563-seed0001-khaase.002.nml', 1 )
%

    % matches .nml format on Feb 27th 2010
    % loads all tags per node (28/06/2010)
   
    % If no path is given, open the console to select a file.
    if nargin < 1
        krk_fname = '';
        [filename, pathname] = uigetfile;
        if filename ~= 0        
            krk_fname = fullfile( pathname, filename );
        end
    end
    
    % Save nodes as struct if not indicated otherwise
    if nargin < 2
        krk_keepNodeAsStruct = 1;
    end
    
    % Load the whole .nml file into Matlab as a huge string
    if ~isempty( krk_fname )
        fid = fopen( krk_fname, 'r' );
        krk_contents = fscanf( fid, '%c' );
        fclose( fid );
        
        
        % Read the part PARAMETERS from krk_contents and split them into
        % several separate cells
        krk_parameters = regexp( krk_contents, '<parameters>.*</parameters>', 'match' );
        krk_thispars = regexp( krk_parameters{1}, '<\w+ [^<]*/>', 'match' );
        krk_parameters_asStruct = struct();
        
        % Write all the separate cells in a struct.
        for krk_parc = 1 : size( krk_thispars, 2 )
            krk_thisparname = strrep( regexp( krk_thispars{ krk_parc }, '<\w+', 'match' ), '<','' );
            krk_thisvalues = regexp( krk_thispars{ krk_parc }, '(\w+)=([^\s/>]+)', 'match' );
            
            for krk_valc = 1 : size( krk_thisvalues, 2 )
                krk_thisNameValues = regexp( krk_thisvalues{ krk_valc }, '=', 'split' );
                krk_parameters_asStruct.( krk_thisparname{1} ).( krk_thisNameValues{1} )...
                    = strrep( krk_thisNameValues{2}, '"', '' );
            end
        end
        
        % Write the paramter-struct into the output struct
        krk_output{1}.parameters = krk_parameters_asStruct;
        
%         krk_parameters = regexp(krk_parameters,'\".*?\"','match');
%         for krk_pc=3:size(krk_parameters{1},2)
%             krk_output{1}.parameters(krk_pc-2) = str2num(krk_parameters{1}{krk_pc}(2:end-1));
%         end
%         krk_output{1}.parameters_dataName = krk_parameters{1}{2};
%         krk_output{1}.parameters_magnification = str2num(krk_parameters{1}{1}(2:end-1));
        

        % Write the part COMMENTS from krk_contents into the output struct
        krk_output{1}.commentsString = regexp( krk_contents, '<comments>.*</comments>', 'match' );
        
        
        % Read the part THING 'ID=X' from krk_contents
        krk_things = regexp( krk_contents, '<thing id.*?</thing>', 'match' );
        
        % Write all nodes and edges into the output struct 
        for krk_tc = 1 : size( krk_things, 2 )
            
            % Write each node individually into a cell and fill up the
            % related space in the output struct with zeros
            %krk_output{krk_tc} = struct('nodes',[],'edges',[]);
            krk_theseNodes = regexp( krk_things{krk_tc}, '<node .*?/>', 'match' );
            krk_output{krk_tc}.nodes = repmat( 0, [size( krk_theseNodes, 2 ), 5] );
            
            % Write now each node successively into the output struct
            for krk_nc = 1 : size( krk_theseNodes, 2 )
                krk_thisNode = regexp( krk_theseNodes{krk_nc}, '\".+?\"', 'match');
                krk_output{krk_tc}.nodes( krk_nc, : ) = ...
                   [str2double( krk_thisNode{1}(2 : end - 1) ),...
                    str2double( krk_thisNode{2}(2 : end - 1) ), ...
                    str2double( krk_thisNode{3}(2 : end - 1) ), ...
                    str2double( krk_thisNode{4}(2 : end - 1) ), ...
                    str2double( krk_thisNode{5}(2 : end - 1) )];
                
                krk_thisNodeTags = regexprep( regexp( krk_theseNodes{krk_nc}, '[ =]', 'split' ),...
                    '["</>]', '' );
                
                % Optionally, write nodes also as struct into the output struct
                if krk_keepNodeAsStruct > 0
                    krk_output{krk_tc}.nodesAsStruct{krk_nc} = struct( krk_thisNodeTags{2 : end} );
                end
                
                % Write only the numerical values of the nodes into a separate part of the
                % output struct
                krk_thisNode_numeric = str2double( krk_thisNodeTags );
                krk_output{krk_tc}.nodesNumDataAll( krk_nc, 1 : length( krk_thisNode_numeric(3 : 2 : end) ) ) = ...
                    krk_thisNode_numeric( 3 : 2 : end );
                
                % In the first run, allocate the matrix for the numerical values with zeros
                if krk_nc == 1
                    krk_output{krk_tc}.nodesNumDataAll( size( krk_theseNodes, 2 ), 1) = 0;
                end
            end
            
            % Create an array with the property: a[node ID] = number of loop cycles,
            % equivalent to a[n] = max(a) - n + 2 for all but a[1] = 1
            krk_nodeIDconversion = zeros( size( krk_theseNodes, 2 ), 1 );
            krk_nodeIDconversion_all = zeros( size( krk_theseNodes, 2 ), 2 );
            
            for krk_nc = 1 : size( krk_theseNodes, 2 )
                krk_nodeIDconversion( krk_output{krk_tc}.nodes( krk_nc, 1 ) ) = krk_nc;
                krk_nodeIDconversion_all( krk_output{krk_tc}.nodes( krk_nc,1 ), 1 : 2) = [krk_tc, krk_nc];
            end
            
            % Write each edge individually into a cell and fill up the
            % related space in the output struct with zeros
            krk_theseEdges = regexp( krk_things{krk_tc}, '<edge .*?/>', 'match' );
            krk_output{krk_tc}.edges = zeros([size( krk_theseEdges, 2 ), 2] );
            
            % Write now each edge successively into the output struct
            for krk_nc = 1 : size( krk_theseEdges, 2 )
                krk_thisEdge = regexp( krk_theseEdges{krk_nc}, '\".+?\"', 'match' );
                krk_output{krk_tc}.edges( krk_nc, : ) =  ...
                   [str2double( krk_thisEdge{1}( 2 : end - 1 ) ), ...
                    str2double( krk_thisEdge{2}( 2 : end - 1 ) )];
            end
            
            % "Invert" the values of .edges 
            krk_output{krk_tc}.edges = krk_nodeIDconversion( krk_output{krk_tc}.edges );
            if size(krk_output{krk_tc}.edges, 2) == 1
                krk_output{krk_tc}.edges = krk_output{krk_tc}.edges';
            end
            % Change the columns of the nodes matrix: Delete [1], put [2] to the end
            krk_output{krk_tc}.nodes = krk_output{krk_tc}.nodes( :, [3 : 5, 2] );
            
            krk_thingID = regexp( krk_things{krk_tc}, '<thing id.{6}', 'match' );
            krk_thingID = regexp( krk_thingID{1}, '[0123456789]*', 'match' );
            krk_output{krk_tc}.thingID = str2double( krk_thingID{1} );
            
        end
        
    end
    

    % Read the part BRANCHPOINTS from krk_contents
    krk_output{1}.branchpointsString = regexp( krk_contents, '<branchpoints>.*</branchpoints>', 'match' );
    krk_output{1}.branchpoints = [];        
    
    if ~isempty( krk_output{1}.branchpointsString )
        % Save only the id of the branchpoints, then convert them to double values
        krk_theseBranchpoints = regexprep( regexp( krk_output{1}.branchpointsString{1}, ...
            '<branchpoint id="', 'split' ), '"/>', '' );
        krk_theseBranchpoints{end} = strrep( krk_theseBranchpoints{end}, '</branchpoints>', '' );
        krk_theseBranchpoints_num = str2double( krk_theseBranchpoints );
        
        % Leave out the first cell, then "invert" the id of the branchpoints
        if length( krk_theseBranchpoints_num ) > 1
            krk_validBranchpoints = krk_theseBranchpoints_num( 2 : end );
            krk_validBranchpoints = krk_validBranchpoints( krk_validBranchpoints <= size( krk_nodeIDconversion_all, 1 ) );
            krk_output{1}.branchpoints = krk_nodeIDconversion_all( krk_validBranchpoints, : );        
        end        
    end

end
