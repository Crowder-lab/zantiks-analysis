#!/bin/sh

# regular analyses
./developmental_delay.py
./distance_travelled.R
./light_dark_preference_3wpf.R
./light_dark_preference_6dpf.R
./light_dark_transition.R
./mirror_biting.R
./sleep.R
./social_preference.R
./startle_response.R
./ymaze_15.R
./ymaze_4.R

# heatmaps
./distance_travelled_heatmap.py
./light_dark_preference_3wpf_heatmap.py
./light_dark_preference_6dpf_heatmap.py
./light_dark_transition_heatmap.py
./mirror_biting_heatmap.py
./sleep_heatmap.py
./social_preference_heatmap.py
./startle_response_heatmap.py
./ymaze_15_heatmap.py
./ymaze_4_heatmap.py
