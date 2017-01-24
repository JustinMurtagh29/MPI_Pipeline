function tests = testGitLabCI()
    tests = functiontests(localfunctions);
end

function testAnswerToLifeAndUniverse(testCase)
    answerToLifeTheUniverseAndEverything = 2 * 3 * 7;
    verifyEqual(testCase, answerToLifeTheUniverseAndEverything, 42);
end