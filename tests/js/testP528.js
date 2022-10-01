// Check if we're running in Node.js (might be running in a browser)
const IS_NODE = (typeof process === 'object' && typeof require === 'function');

// We're in the test runner (we'll load these in the 'before' function)
if (IS_NODE) {
    let chai = null;
    let Module = null;
}
else { // We're in the browser...
    // Using 'var' because it needs to be global in order for Emscripten's generated JavaScript code to see it
    var Module = {
        // Run the tests when Emscripten's JS tells us that it's ready to be interacted with
        onRuntimeInitialized: () => { mocha.run(); }
    };
}

describe("Testing the model P.528", () => {

    before(() => {
        if (IS_NODE) {
            // Load the Chai library
            chai = require('chai');

            // Load the Emscripten generated JavaScript. Because it might not be ready right away, return a Promise that will be
            // resolved once the module is ready to be interacted with
            return new Promise((resolve) => {
                Module = require('../../wasm/dist/p528json');
                Module['onRuntimeInitialized'] = () => {
                    resolve();
                }
            }); // new Promise
        } // IS_NODE
    }); // before

    it("Calculate Profile", async () => {

        const request = {
            h_1__user_units: 5,
            h_2__user_units: 10000,
            time: 50,
            Polarization: 0,
            f__mhz: 5000, 
            nbr_of_points: 5
        }

        const requestStr = JSON.stringify(request);

        console.log("request: " + requestStr + "\n");

        const outputJson = Module.ccall('calcP528', 'string', ['string'], [requestStr]);

        console.log("output: " + outputJson + "\n");

        const expectedJson = "{\"A__db\":[126.46085394936013,163.0815203486053,236.28197179677775,268.19292392138044,295.0633534171327,323.82517691745284],\"A_a__db\":0,\"A_fs__db\":[126.425056055803,157.56187253649026,163.55576213633907,167.04435177821983,169.5660332409185,171.5603915860845],\"ModeOfPropagation\":[1,1,3,3,3,3],\"d__km\":[0.0,360.0008519649057,720.0,1080.0,1440.0,1800.0],\"theta_h1__rad\":0,\"warnings\":0}";

        chai.expect(expectedJson).to.equal(outputJson);
    });
});
