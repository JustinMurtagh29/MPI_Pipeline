#!/usr/bin/env bash
echo "> Running tests..."
export PATH=$PATH:/usr/local/matlab/r2017b/bin
matlab -nosplash -nodisplay -r "runAllTests();"
