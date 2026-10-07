# Author: Mazlina Ismail
# Based on ARN86.R

load('/path/to/2023-10-16-ARNVER55-arneo-array-data-noLN.RData')
# select continuous sig scores only and keep patient ID
keep_dat <- scores_dat[,c(469, 470, 2:468)]
# PY <- biopsy
py <- subset(keep_dat, label == 'PY')
py_mean <- aggregate(py[,3:ncol(py)], by=list(PATIENT_ID=py$PATIENT_ID), FUN=mean)
# PR <- prostatectomy
pr <- subset(keep_dat, label == 'PR')
pr_mean <- aggregate(pr[,3:ncol(pr)], by=list(PATIENT_ID=pr$PATIENT_ID), FUN=mean)
