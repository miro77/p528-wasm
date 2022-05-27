mkdir -p ./public/js/
rm -rf ./public/js/*
cp ../wasm/dist/p528json.wasm ./public/js/ 
cp ../wasm/dist/p528json.js ./public/js/
cp ../wasm/dist/p528json.worker.js ./public/js/
