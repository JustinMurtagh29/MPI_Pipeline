#!/usr/bin/env bash
echo "> Connecting to license server..."
ssh -fN -L 27000:172.16.1.10:27000 -L 28000:172.16.1.10:28000 gaba.opt.rzg.mpg.de
