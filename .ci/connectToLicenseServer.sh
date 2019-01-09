#!/usr/bin/env bash
echo "> Connecting to license server..."
ssh -4 -fN -L 27000:192.109.27.61:27000 -L 28000:192.109.27.61:28000 gaba.opt.rzg.mpg.de
