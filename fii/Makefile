#########################################################################################
# Copyright 2020 Konstantinos Zagganas and Thanasis Vergoulis
# for the Information Management Systems Institute(IMSI) - "Athena" Research Center
# 
# This file is part of diwis/fii-sli.
#
# Foobar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# diwis/fii-sli is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
#########################################################################################
CC=g++

CFLAGS=-std=c++14 -O3 -static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

fii: bufet-fii.cpp
	$(CC) bufet-fii.cpp -o bufet-fii $(CFLAGS)