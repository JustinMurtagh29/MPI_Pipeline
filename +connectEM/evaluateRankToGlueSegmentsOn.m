function yy = evaluateRankToGlueSegmentsOn(ygrid5, segmentMeta,excludedSegmentIdx,gridAgglo_05, graph, agnostic)
    function y = nonull(x)
        y = unique(x(x~=0));
    end
    tempbbb = cellfun(@(x)feval(@(y)y(excludedSegmentIdx(y)), nonull(x)),ygrid5{564}.axon1.ids, 'uni', false);
    tempc = cellfun(@(x)cell2mat(gridAgglo_05{564}.axonsFinal(x)), ygrid5{564}.axon1.foundAgglomerates_col, 'uni', 0);
    tempddd = cellfun(@(x,y)setdiff(nonull(x), y),ygrid5{564}.axon1.ids,tempc, 'uni', 0);

    tempaaa = cellfun(@(x)feval(@(y)y(segmentMeta.axonProb(y)<0.5), nonull(x)),ygrid5{564}.axon1.ids, 'uni', 0);

    for idx = 1 : 10
        idx

        if agnostic
            tempfff{idx} = setdiff(tempaaa{idx}, [tempbbb{idx}]);
        else
            tempfff{idx} = setdiff(tempddd{idx}, [tempbbb{idx}; tempaaa{idx}]);
        end
        yy.VV{idx} = [];
        yy.VVinfo{idx} = [];
        yy.VVrank{idx} = [];
        yy.VVrankInfo{idx} = {};
        for idx2 = 1 : length(tempfff{idx})
            connectedAgglos = cellfun(@(x)[-1; graph.neighProb{tempfff{idx}(idx2)}(ismember(graph.neighbours{tempfff{idx}(idx2)},x))], gridAgglo_05{564}.axonsFinal(ygrid5{564}.axon1.foundAgglomerates_col{idx}), 'uni', 0);
            [V, I] = max(cellfun(@max, connectedAgglos));
            if V >= 0
                thisagglo = gridAgglo_05{564}.axonsFinal{ygrid5{564}.axon1.foundAgglomerates_col{idx}(I)};
                tann3 = setdiff(cell2mat(graph.neighbours(thisagglo)'), [thisagglo; tempfff{idx}(idx2)]);
                tann3(excludedSegmentIdx(tann3)) = [];
                if ~agnostic
                    tann3(segmentMeta.axonProb(tann3)<0.5) = [];
                end
                clear tann
                for idx3 = 1 : length(tann3)
                    tann(idx3) = max(cell2mat(cellfun(@(x, y) max(y(x==tann3(idx3))), graph.neighbours(thisagglo),graph.neighProb(thisagglo), 'uni', 0)));
                end
                yy.VV{idx}(end + 1) = V;
                yy.VVinfo{idx}(end + 1) = tempfff{idx}(idx2);
                yy.VVrank{idx}(end + 1) = sum(tann > V);
                yy.VVrankInfo{idx}{end + 1} = tann3(tann > V);
            end
        end
    end
end
