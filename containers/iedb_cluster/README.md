# Info

Available prediction thresholds are:
------------------------------------
0 to 100 (default=80)

Available values for unique option are:
---------------------------------------
'on' or 'off' (default=on)

Available commands:
-------------------
python predict_cluster.py [input-file]
python predict_cluster.py [options] [input-file]

Example:
--------
python predict_cluster.py example/test.txt
python predict_cluster.py --threshold=70 --unique=off example/test.txt
