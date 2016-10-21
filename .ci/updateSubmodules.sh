#!/usr/bin/env bash
echo "> Updating submodules..."
echo -e "[url \"git@gitlab.mpcdf.mpg.de:\"]\n insteadOf = \"https://gitlab.mpcdf.mpg.de/\"" > ~/.gitconfig
git submodule update --init --recursive
