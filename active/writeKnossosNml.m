function writeKnossosNml( ksw_fname, ksw_skeleton )

% WRITEKNOSSOSNML: Write Knossos skeletons in Matlab as a .nml file
%   
%   The function has the following arguments:
%       KSW_FNAME: Give the full path of the file in which you want to
%           write the data, e.g. 'E:\knossos_skeleton.nml'
%       KSW_SKELETON: Give the name of the cell array containing the
%           Knossos skeleton(s), e.g. tracing in the Matlab
%           workspace.
%       
%   => writeKnossosNml( 'E:\knossos_skeleton.nml', tracing )
%

    % Open the file, if it already exists, overwrite the contents.
    fid = fopen( ksw_fname, 'w' );
    
    % The .nml format is an ASCII file, thus human readable.
    fprintf( fid, '<?xml version=\"1.0\"?>\n' );
    fprintf( fid, '<things>\n' );

    % If ksw_skeleton{1}.parameters is a structure array, continue.
    if isfield( ksw_skeleton{1}, 'parameters')
        
        % Read out the necessary information.
        ksw_experimentName = ksw_skeleton{1}.parameters.experiment.name;
        ksw_scale = ksw_skeleton{1}.parameters.scale;
        ksw_offset = ksw_skeleton{1}.parameters.offset;
        ksw_time = ksw_skeleton{1}.parameters.time.ms;
        ksw_activeNode = ksw_skeleton{1}.parameters.activeNode.id;
        ksw_editPosition = ksw_skeleton{1}.parameters.editPosition;
        
        % Write the parameters into the file.
        fprintf( fid, '\t<parameters>\n\t\t<experiment name=\"%s\"/>\n\t\t<scale x=\"%s\" y=\"%s\" z=\"%s\"/>\n\t\t<offset x=\"%s\" y=\"%s\" z=\"%s\"/>', ...
        ksw_experimentName, ksw_scale.x, ksw_scale.y, ksw_scale.z, ...
        ksw_offset.x, ksw_offset.y, ksw_offset.z );
    
        fprintf( fid, '\n\t\t<time ms=\"%s\"/>\n\t\t<activeNode id=\"%s\"/>\n\t\t<editPosition x=\"%s\" y=\"%s\" z=\"%s\"/>\n\t</parameters>\n', ...
            ksw_time, ksw_activeNode, ksw_editPosition.x, ksw_editPosition.y, ksw_editPosition.z );
        clear( 'ksw_experimentName', 'ksw_scale', 'ksw_offset', 'ksw_time', 'ksw_activeNode', 'ksw_editPosition');
    end
    
    % Necessary if multiple skeletons exist.
    kl_totalNodeC = 0;
    
    % Determine the number of different skeletons and go over each one.
    for kl_thingC = 1 : size( ksw_skeleton, 2 )
        
        % Start with writing the nodes into the file.
        fprintf( fid, '\t<thing id=\"%d\">\n\t\t<nodes>\n', kl_thingC );
        
        if isfield( ksw_skeleton{kl_thingC}, 'nodesNumDataAll' ) 
            
            ksw_nodeList = ksw_skeleton{kl_thingC}.nodesNumDataAll;
            
            for ksw_c = 1 : size( ksw_nodeList, 1 )
                fprintf( fid, '\t\t\t<node id=\"%d\" radius=\"%.6f\" x=\"%d\" y=\"%d\" z=\"%d\" inVp=\"%d\" inMag=\"%d\" time=\"%d\"/>\n', ...
                    ksw_c + kl_totalNodeC, ksw_nodeList( ksw_c, 2 ), ksw_nodeList( ksw_c, 3 ), ...
                    ksw_nodeList( ksw_c, 4 ), ksw_nodeList( ksw_c, 5 ), ksw_nodeList( ksw_c, 6 ), ...
                    ksw_nodeList( ksw_c, 7 ), ksw_nodeList( ksw_c, 8 ) );
            end
            
        else
          
            ksw_nodeList = ksw_skeleton{kl_thingC}.nodes;
        
            % Write each node successively into the file.
            for ksw_c = 1 : size( ksw_nodeList, 1 )
                fprintf( fid, '\t\t\t<node id=\"%d\" radius=\"%.6f\" x=\"%d\" y=\"%d\" z=\"%d\"/>\n', ...
                    ksw_c + kl_totalNodeC, ksw_nodeList( ksw_c, 4 ), ksw_nodeList( ksw_c, 1 ), ...
                    ksw_nodeList( ksw_c, 2 ), ksw_nodeList( ksw_c, 3 ) );
            end
        end
        
        fprintf( fid, '\t\t</nodes>' );
        
        % Start with writing the edges into the file. If the matrix does
        % not exist or is empty, simply write 'edges' into the file.
        if ~isfield( ksw_skeleton{kl_thingC}, 'edges' ) || isempty( ksw_skeleton{kl_thingC}.edges )
            fprintf( fid, '\n\t\t<edges/>' );
        
        % If edges is a structure array, write the edges successively into
        % the file.
        else
            fprintf( fid, '\n\t\t<edges>\n' );
            
            for ksw_ec = 1 : size( ksw_skeleton{kl_thingC}.edges, 1 )
                fprintf( fid, '\t\t\t<edge source=\"%d\" target=\"%d\"/>\n', ...
                    ksw_skeleton{kl_thingC}.edges( ksw_ec, 1 ) + kl_totalNodeC, ...
                    ksw_skeleton{kl_thingC}.edges( ksw_ec, 2 ) + kl_totalNodeC );
            end
            
            fprintf( fid, '\n\t\t</edges>' );
        end

        % Change kl_totalNodeC to let the node id start at the right number
        % if multiple skeletons exist.
        kl_totalNodeC = kl_totalNodeC + size( ksw_nodeList, 1 );
        fprintf( fid, '\n\t</thing>\n' );
    end
    
    % Write the comments into the file (new version of saving). The node count
    % in the comments and the node id will match even with multiple tracings.
    if iscell( ksw_skeleton{1}.commentsString )
        
        kl_totalNodeC = 0;
        fprintf( fid, '\t<comments>\n' );
        
%         % Go over all the tracings.
%         for kl_thingC = 1 : size( ksw_skeleton, 2 )
%             
%             % Extract each single comment and write it into a cell array.
%             ksw_comments = ksw_skeleton{1}.commentsString{1}( 17 : ( end - 15) );
%             ksw_comments = regexp( ksw_comments, '<\w+ [^<]*/>', 'match' );
%             
%             for ksw_c = 1 : size( ksw_comments , 2)
%                 
%                 % Extract the node id from the comments string.
%                 ksw_nodeId = regexp( ksw_comments{ksw_c}, '"\d+"', 'match' );
%                 ksw_nodeId = str2double( ksw_nodeId{1}( 2 : ( end - 1 ) ) );
%                 
%                 % Convert the node id into the current one.
%                 if ksw_nodeId ~= 1    
%                     ksw_nodeId = size( ksw_skeleton{kl_thingC}.nodes, 1 ) - ( ksw_nodeId - 2 );
%                     ksw_nodeId = ksw_nodeId + kl_totalNodeC;                    
%                 end
%                 
%                 if ksw_nodeId == 1 && kl_totalNodeC ~= 0
%                    ksw_nodeId = kl_totalNodeC + 1; 
%                 end
%                 
%                 % Extract the comment from whole comments string.
%                 ksw_comment = regexp( ksw_comments{ksw_c}, '"[ a-zï¿½ï¿½ï¿½ï¿½A-Zï¿½ï¿½ï¿½]+"', 'match' );
%                 
%                 % Write the comment into the file.
%                 if ~isempty(ksw_comment)
%                     fprintf( fid, '\t\t<comment node=\"%d\" content=%s/>\n',...
%                         ksw_nodeId, ksw_comment{1} );
%                 end
%                 
%             end
%             
%             % Increase kl_totalNodeC analogously.
%             kl_totalNodeC = kl_totalNodeC + size( ksw_skeleton{kl_thingC}.nodes, 1 );
%         end 
%         
%         fprintf( fid, '\t</comments>\n' );
%     end
%     
%     % Write the comments into the file (old version of saving)
%     if isfield( ksw_skeleton{1}, 'comments' )
%         fprintf( fid, '\t<comments>\n' );
%         
%         for kl_cc = 1 : size( ksw_skeleton{1}.comments, 1 )
%             fprintf( fid, '\t\t<comment node=\"%d\" content=\"%s\"/>\n', ...
%                  ksw_skeleton{1}.comments{ kl_cc, 1 }, ksw_skeleton{1}.comments{ kl_cc, 2 } );
%         end
%         
%         fprintf( fid, '\t</comments>\n' );
%     end
    fprintf( fid, '\t</comments>\n' );
    % Write the last line, then close the file. 
    fprintf( fid, '</things>\n' );
    fclose( fid );
    
    % Print on the screen.
    fprintf( 'Done writing!\n' );
    
end
