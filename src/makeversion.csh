#!/bin/csh

set ts=(`date`)

echo 'static char Version_Label[]={"3.0.5"}; ' >! ../include/mt_version.h

echo 'static char Version_Date[]={"\' >> ../include/mt_version.h
echo "${ts}\" >> ../include/mt_version.h
echo '" };' >> ../include/mt_version.h
