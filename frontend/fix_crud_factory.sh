#!/bin/bash
cd /home/user/NGSmodule/frontend

# Fix all occurrences where _state parameter is used but state is referenced in body
# Use sed to replace the pattern, being careful to only change the parameter name
sed -i '312s/set((_state)/set((state)/' src/store/crud.factory.ts
sed -i '350s/set((_state)/set((state)/' src/store/crud.factory.ts  
sed -i '412s/set((_state)/set((state)/' src/store/crud.factory.ts

# Also need to fix the get parameter
sed -i '198s/, _get)/, get)/' src/store/crud.factory.ts
sed -i '449s/get: () => BaseState & BaseActions & CustomState & CustomActions/_get: () => BaseState & BaseActions & CustomState & CustomActions/' src/store/crud.factory.ts
sed -i '456s/get as () => State & Actions/_get as () => State & Actions/' src/store/crud.factory.ts
sed -i '483s/get as () => State & Actions/_get as () => State & Actions/' src/store/crud.factory.ts

echo "Fixed crud.factory state/get parameters"
