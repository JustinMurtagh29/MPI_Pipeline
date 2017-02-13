function rewriteWithoutOverlaps( inFile, outFile, bbS, bbB)

    seg = extractBboxSmall(inFile, bbS, bbB);
    numEl = calcNumberSegments(seg);
    Util.save(outFile, seg, numEl);

end

