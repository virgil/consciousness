#!/bin/sh

hadoop jar /usr/local/hadoop-0.19.0/contrib/streaming/hadoop-0.19.0-streaming.jar -input /todo -output /todone -mapper /consciousness/consciousness -reducer /consciousness/MR-reducer.py
