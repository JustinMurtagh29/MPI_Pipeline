function tests = testGitLabCI()
    tests = functiontests(localfunctions);
end

function testAnswerToLifeAndUniverse(testCase)
    answerToLifeAndUniverse = 2 * 3 * 7;
    verifyEqual(testCase, answerToLifeAndUniverse, 42);
end