#!/bin/bash
cd /home/user/NGSmodule/frontend

# Fix unused variables that sed didn't catch properly
# ScatterPlot
sed -i 's/  showRegression,$/  _showRegression,/' src/components/charts/ScatterPlot.tsx

# ConfirmDialog
sed -i 's/  type,$/  _type,/' src/components/common/ConfirmDialog.tsx

# PageSkeleton
sed -i 's/  columns,$/  _columns,/' src/components/common/Skeleton/PageSkeleton.tsx

# BatchImportModal
sed -i 's/(_, record)/(_, _record)/' src/components/samples/BatchImportModal.tsx

# TaskProgressCard - Title
sed -i 's/const { Title }/const { Title: _Title }/' src/components/tasks/TaskProgressCard.tsx

# FileList - info parameter
sed -i 's/info: FileInfo/  _info: FileInfo/' src/pages/files/FileList.tsx

# FileList - handleUpload function (it's actually not used, remove it)
sed -i '/_handleUpload/d' src/pages/files/FileList.tsx

# ResultList - error variable
sed -i '74s/error/  _error/' src/pages/results/ResultList.tsx

# project.service.ts - PaginatedResponse
sed -i '/^import type { PaginatedResponse }/d' src/services/project.service.ts

# stats.service.ts - taskStats
sed -i 's/const taskStats/const _taskStats/' src/services/stats.service.ts

# websocket.service.ts - url
sed -i 's/const url =/const _url =/' src/services/websocket.service.ts

# crud.factory.ts - state parameter
sed -i 's/(state) => (/(_state) => (/' src/store/crud.factory.ts

echo "Fixed unused variables"
