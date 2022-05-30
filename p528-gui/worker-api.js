import Worker from 'worker-loader!./worker'

import * as Comlink from 'comlink'

const worker = Comlink.wrap(new Worker('./worker'), { type: 'module' })

export const calc = worker.calc