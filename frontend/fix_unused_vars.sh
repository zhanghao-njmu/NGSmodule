#!/bin/bash

# Fix unused variables by prefixing with underscore

# ScatterPlot.tsx
sed -i 's/  showRegression,/  _showRegression,/' src/components/charts/ScatterPlot.tsx

# ConfirmDialog.tsx
sed -i 's/  type,/  _type,/' src/components/common/ConfirmDialog.tsx

# PageSkeleton.tsx
sed -i 's/  columns,/  _columns,/' src/components/common/Skeleton/PageSkeleton.tsx

# BatchImportModal.tsx
sed -i 's/, record)/, _record)/' src/components/samples/BatchImportModal.tsx

# TaskProgressCard.tsx
sed -i 's/const { Title } =/const { Title: _Title } =/' src/components/tasks/TaskProgressCard.tsx

# Dashboard.tsx - remove unused useState
sed -i '/^import { useState }/d' src/pages/dashboard/Dashboard.tsx

# FileList.tsx
sed -i 's/info:/  _info:/' src/pages/files/FileList.tsx
sed -i 's/const handleUpload =/const _handleUpload =/' src/pages/files/FileList.tsx

# ProjectFormModal.tsx - remove unused message
sed -i 's/, message } from/@ant-design\/icons/' src/pages/projects/components/ProjectFormModal.tsx

# ResultList.tsx - will handle separately (duplicate error variable)

# ai.service.ts
sed -i 's/callback:/  _callback:/' src/services/ai.service.ts

# stats.service.ts
sed -i 's/const taskStats =/const _taskStats =/' src/services/stats.service.ts

# websocket.service.ts
sed -i 's/const url =/const _url =/' src/services/websocket.service.ts

# authStore.ts
sed -i 's/const user =/const _user =/' src/store/authStore.ts

# crud.factory.ts
sed -i 's/(set, get)/(set, _get)/' src/store/crud.factory.ts
sed -i "s/set\((state)\)/set\((_state)\)/" src/store/crud.factory.ts

echo "Unused variables fixed"
