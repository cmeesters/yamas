#! /bin/bash

# Copyright (C) 2010, Christian Meesters (meesters@imbie.uni-bonn.de)
# This code is part of the bimp distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# cleanup the doc directory after a LaTeX run

# $Rev: 785 $
# $LastChangedDate: 2011-02-13 19:30:52 +0100 (So, 13 Feb 2011) $ # date of last revision
# $Author: meesters $

\rm -f ./doc/*.blg
\rm -f ./doc/*.out
\rm -f ./doc/*.aux
\rm -f ./doc/*~
\rm -f ./doc/*.log
\rm -f ./doc/*.toc
\rm -f ./doc/*.bbl
\rm -f ./doc/*.backup
\rm -f ./paper/*.blg
\rm -f ./paper/*.out
\rm -f ./paper/*.aux
\rm -f ./paper/*~
\rm -f ./paper/*.log
\rm -f ./paper/*.toc
\rm -f ./paper/*.bbl
\rm -f ./paper/*.backup

