classdef BoostedEnsemble
    %ENSEMBLE Wrapper class for MATLABs ensembles (mainly AdaBoostM1 and
    % LogitBoost.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        ens             %matlab ensemble class
        options         %training options
        predVar = []    %variables required for prediction
    end
    
    methods
        function [y, scores] = predict(obj, X, learners)
            %PREDICT Default prediction function for ens if predVar is
            %empty. Otherwise uses pred var for faster prediction which
            %produces exactly the same result.
            % INPUT X: [NxM] float feature matrix.
            %       learners: (Optional) See ens.predict.
            % OUTPUT y: [Nx1] logical. Prediction for each row in X.
            %        scores: [Nx1] float containing the scores for each row
            %           in y. High scores correspond to higher likelihood
            %           of the example to be in class true.
            
            %faster prediction for decision stumps
            if ~isempty(obj.predVar)
                if ~exist('learners','var') || isempty(learners)
                    learners = 1:length(obj.predVar.W);
                end
                scores = zeros(size(X,1),1);
                for i = learners
                    indices = X(:,obj.predVar.cutVars(i)) ...
                        < obj.predVar.cutPoints(i);
                    scores = scores + obj.predVar.W(i).* ...
                        ((indices.*obj.predVar.leftNodePred(i) + ...
                        ~indices.*obj.predVar.rightNodePred(i)));
                end
                y = scores > 0;
            else
                if ~exist('learners','var') || isempty(learners)
                    [y, scores] = predict(obj.ens, X);
                else
                    [y, scores] = predict(obj.ens, X, ...
                        'learners', learners);
                end
                scores = scores(:,1);
            end
                
        end
        
        function obj = calculatePredVar(obj)
            % CALCULATEPREDVAR Extract the variables used for 
            % prediction from an ensemble of decision stumps for a
            % binary classifier.
            % NOTE Only use for binary classification and decision stumps.
            %      (Only enforces via some sanity checks below).
            
            if length(obj.ens.ClassNames) ~= 2
                error(['Currently only implemented for binary', ...
                    'classification.']);
            end
            
            if strcmp(obj.options.method, 'LogitBoost')
                if obj.ens.Trained{1}.CompactRegressionLearner.NumNodes ~= 3
                    error(['Currently only implemented for trees with', ...
                        'one split (default for fitensemble).']);
                end
                T = obj.ens.NumTrained;
                W = obj.ens.TrainedWeights;
                cutVars = zeros(T,1);
                cutPoints = zeros(T,1);
                leftNodePred = zeros(T,1);
                rightNodePred = zeros(T,1);
                for t = 1:T
                    weakLearner = obj.ens.Trained{t}.CompactRegressionLearner;
                    cutPoints(t) = weakLearner.CutPoint(1);
                    cutVars(t) = str2double(weakLearner.CutVar{1}(2:end));
                    leftNodePred(t) = -weakLearner.NodeMean(2);
                    rightNodePred(t) = -weakLearner.NodeMean(3);
                end
                obj.predVar.W = W;
                obj.predVar.cutVars = cutVars;
                obj.predVar.cutPoints = cutPoints;
                obj.predVar.leftNodePred = leftNodePred;
                obj.predVar.rightNodePred = rightNodePred;
            elseif strcmp(obj.options.method, 'AdaBoostM1')
                if obj.ens.Trained{1}.NumNodes ~= 3
                    error(['Currently only implemented for trees with', ...
                        'one split (default for fitensemble).']);
                end
                T = obj.ens.NumTrained;
                W = obj.ens.TrainedWeights;
                cutVars = zeros(T,1);
                cutPoints = zeros(T,1);
                leftNodePred = zeros(T,1);
                rightNodePred = zeros(T,1);
                for t = 1:T
                    weakLearner = obj.ens.Trained{t};
                    cutPoints(t) = weakLearner.CutPoint(1);
                    cutVars(t) = str2double(weakLearner.CutVar{1}(2:end));
                    if strcmp(weakLearner.NodeClass{2},'true')
                        leftNodePred(t) = 1;
                    else
                        leftNodePred(t) = -1;
                    end
                    if strcmp(weakLearner.NodeClass{3},'true')
                        rightNodePred(t) = 1;
                    else
                        rightNodePred(t) = -1;
                    end
                end
                obj.predVar.W = W;
                obj.predVar.cutVars = cutVars;
                obj.predVar.cutPoints = cutPoints;
                obj.predVar.leftNodePred = leftNodePred;
                obj.predVar.rightNodePred = rightNodePred;
            else
                error(['Currently not implemented for ', ...
                    'classification method %s.'], obj.ens.Method);
            end
        end
        
        function obj = compact(obj)
            %Call the compact function from MATLABs ensemble class on ens.
            obj.ens = compact(obj.ens);
        end
    end
    
    methods (Static)
        function obj = train(X, y, varargin)
            %Wrapper for fitensemble with some default options.
            % INPUT X: [NxM] float containing the feature matrix. Rows
            %           correspond to observations and columns to
            %           variables.
            %       y: [Nx1] logical containing the labels for the rows in
            %           X.
            %       varargin: Name value pairs for the fitensemble options.
            
            t = datestr(clock,30);
            options = struct;
            options.method = 'LogitBoost';
            options.nlearn = 1500;
            options.learners = 'tree';
            options.LearnRate = 0.1;
            options.nprint = 100;
            options.prior = 'empirical';
            options.classname = [true,false];
            options.crossval = 'off';
            options.cost = [0 100;1 0];
            options.type = 'classification';
            options = Util.modifyStruct(options, varargin{:});
            nmPairs = Util.struct2nmcell(options);
            
            fprintf('[%s] Starting classifier training.\n', datestr(now));
            obj = SynEM.Classifier.BoostedEnsemble();
            obj.ens = fitensemble(X,y,options.method,...
                            options.nlearn,options.learners, ...
                            nmPairs{7:end});
            fprintf('[%s] Finished classifier training.\n', datestr(now));
            options.start = t;
            obj.options = options;
            try %#ok<TRYNC>
                obj = obj.calculatePredVar();
            end
        end
    end
    
end

