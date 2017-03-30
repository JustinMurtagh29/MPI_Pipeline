function comment = makeComment(comment, exclusions, newid)
for fn = fieldnames(exclusions)'
    if exclusions.(fn{:})(newid)
        comment = [comment '_' fn{:}];
    end
end
