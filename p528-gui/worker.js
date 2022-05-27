import * as Comlink from 'comlink'

self.importScripts('p528json.js')

const obj = {
    calc(intputJson) {

        var resultString = self.Module.ccall(
            "calcP528",
            "string",
            ["string"],
            [intputJson]
        );
        return resultString;
    }
};

Comlink.expose(obj);