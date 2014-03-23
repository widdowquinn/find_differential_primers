************* Module stitch_six_frame_stops
C:  1, 0: Missing module docstring (missing-docstring)
C:154, 4: Invalid constant name "options" (invalid-name)
C:154,13: Invalid constant name "args" (invalid-name)
C:158, 4: Invalid constant name "logger" (invalid-name)
C:160, 4: Invalid constant name "err_handler" (invalid-name)
C:161, 4: Invalid constant name "err_formatter" (invalid-name)
C:176, 8: Invalid constant name "inhandle" (invalid-name)
C:179, 8: Invalid constant name "inhandle" (invalid-name)
C:180, 4: Invalid constant name "data" (invalid-name)
C:184, 4: Invalid constant name "data" (invalid-name)
C:190, 4: Invalid constant name "stitchedseq" (invalid-name)
C:195, 8: Invalid constant name "outhandle" (invalid-name)
C:197, 8: Invalid constant name "outhandle" (invalid-name)


Report
======
79 statements analysed.

Messages by category
--------------------

+-----------+-------+---------+-----------+
|type       |number |previous |difference |
+===========+=======+=========+===========+
|convention |13     |13       |=          |
+-----------+-------+---------+-----------+
|refactor   |0      |0        |=          |
+-----------+-------+---------+-----------+
|warning    |0      |0        |=          |
+-----------+-------+---------+-----------+
|error      |0      |0        |=          |
+-----------+-------+---------+-----------+



Messages
--------

+------------------+------------+
|message id        |occurrences |
+==================+============+
|invalid-name      |12          |
+------------------+------------+
|missing-docstring |1           |
+------------------+------------+



Global evaluation
-----------------
Your code has been rated at 8.35/10 (previous run: 8.35/10, +0.00)

Raw metrics
-----------

+----------+-------+------+---------+-----------+
|type      |number |%     |previous |difference |
+==========+=======+======+=========+===========+
|code      |103    |54.50 |103      |=          |
+----------+-------+------+---------+-----------+
|docstring |10     |5.29  |10       |=          |
+----------+-------+------+---------+-----------+
|comment   |64     |33.86 |64       |=          |
+----------+-------+------+---------+-----------+
|empty     |12     |6.35  |12       |=          |
+----------+-------+------+---------+-----------+



Statistics by type
------------------

+---------+-------+-----------+-----------+------------+---------+
|type     |number |old number |difference |%documented |%badname |
+=========+=======+===========+===========+============+=========+
|module   |1      |1          |=          |0.00        |0.00     |
+---------+-------+-----------+-----------+------------+---------+
|class    |0      |0          |=          |0           |0        |
+---------+-------+-----------+-----------+------------+---------+
|method   |0      |0          |=          |0           |0        |
+---------+-------+-----------+-----------+------------+---------+
|function |3      |3          |=          |100.00      |0.00     |
+---------+-------+-----------+-----------+------------+---------+



Duplication
-----------

+-------------------------+------+---------+-----------+
|                         |now   |previous |difference |
+=========================+======+=========+===========+
|nb duplicated lines      |0     |0        |=          |
+-------------------------+------+---------+-----------+
|percent duplicated lines |0.000 |0.000    |=          |
+-------------------------+------+---------+-----------+



External dependencies
---------------------
::

    Bio 
      \-Seq 
      | \-Seq (stitch_six_frame_stops)
      \-SeqIO (stitch_six_frame_stops)
      \-SeqRecord 
        \-SeqRecord (stitch_six_frame_stops)



