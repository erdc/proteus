#!/bin/bash +x

/lore/zhanga/core-sim/build/test/convert $1.smd $1.sms $1.smb
/lore/zhanga/core-sim/build/test/mdlConvert $1.smd $1.dmg
