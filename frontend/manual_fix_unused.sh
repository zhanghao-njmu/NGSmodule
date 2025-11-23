#!/bin/bash
cd /home/user/NGSmodule/frontend

# ScatterPlot - line 30, destructure showRegression but prefix it  
perl -i -pe 's/^(\s+)showRegression,/$1_showRegression,/' src/components/charts/ScatterPlot.tsx

# ConfirmDialog - line 30, destructure type but prefix it
perl -i -pe 's/^(\s+)type,/$1_type,/' src/components/common/ConfirmDialog.tsx

# PageSkeleton - line 19, destructure columns but prefix it
perl -i -pe 's/^(\s+)columns,/$1_columns,/' src/components/common/Skeleton/PageSkeleton.tsx

# BatchImportModal - line 284, parameter record -> _record
sed -i 's/(_, record)/(_, _record)/' src/components/samples/BatchImportModal.tsx

# TaskProgressCard - line 27, Title -> _Title
sed -i 's/const { Title }/const { Title: _Title }/' src/components/tasks/TaskProgressCard.tsx

# FileList - line 76, info parameter
sed -i 's/onChange(info)/onChange(_info)/' src/pages/files/FileList.tsx

# stats.service.ts - line 103, taskStats
sed -i 's/const taskStats =/const _taskStats =/' src/services/stats.service.ts

# websocket.service.ts - line 15, url
sed -i 's/const url =/const _url =/' src/services/websocket.service.ts

echo "Manual fixes applied"
