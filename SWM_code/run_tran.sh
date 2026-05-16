#!/bin/bash

cd main

sh verintp.pres2sigma.sh
sh horzintp.sh
sh ps_horzintp.sh
sh make_basicstate.sh

cd ..
sh tran_forcing.sh
