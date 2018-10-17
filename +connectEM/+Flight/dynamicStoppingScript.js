/* dynamicStoppingScript.js
 *
 * This JavaScript snippet uses to webKNOSSOS API to intersect the nodes of a
 * flight mode tracing against an invisible segmentation while they are added.
 * Once there is enough evidence for a true overlap with an agglomerate
 * different from the start location, the script notifies the user, so that she
 * / he can jump to the next one or take a break.
 *
 * Added to this repository from
 *   https://gist.github.com/heikowissler/1e715e790ac52a11bd24cc8e050bd1fd
 *
 * Written by
 *   Kevin M. Boergens <kevin.boergens@brain.mpg.de>
 *   Alessandro Motta <alessandro.motta@brain.mpg.de>
 */
webknossos.apiReady(2).then(async (api) => {
    console.log("starting jumpy queries");
    window.collection = [];
    window.starter = [];
    window.taskCounter = 0;
    window.done = false;
    api.utils.showToast("info", 'Press "7" to automatically place a comment for unsolved tracing (means: mistracing, glia, unsure) and "9" to finish and get next task', 6000);
    async function createNodeOverwrite(store, call, action) {
        const taskCounterHere = window.taskCounter;
        call(action);
        if (window.collection[taskCounterHere] === undefined) {
            window.collection[taskCounterHere] = [];
        }
        const numberNodes = api.tracing.getAllNodes().length;
        const pos = action.position;
        const posx = pos[0];
        const posy = pos[1];
        const posz = pos[2];
        for (var idx = -1; idx < 2; idx++) {
            for (var idy = -1; idy < 2; idy++) {
                for (var idz = -1; idz < 2; idz++) {
                    const segmentId = await api.data.getDataValue("segmentation", [posx + idx, posy + idy, posz + idz]);
                    if (taskCounterHere < window.taskCounter) {
                        return;
                    }
                    if (segmentId > 0) {
                        if (window.collection[taskCounterHere][segmentId] === undefined) {
                            window.collection[taskCounterHere][segmentId] = 1;
                        } else {
                            window.collection[taskCounterHere][segmentId] += 1;
                            if (window.collection[taskCounterHere][segmentId] >= 54 && (numberNodes > 5 || segmentId !== window.starter[taskCounterHere])) {
                                // task finished
                                api.utils.showToast("success", 'You have finished this task, thanks a lot! You can press "9" now.', 8000);
                                // var layers = api.data.getConfiguration("layers");
                                // layers.color.color = [0, 255, 0];
                                // api.data.setConfiguration("layers", layers);
                                window.done = true;
                            }
                        }
                    }
                }
            }
        }
    }
    async function wKReadyOverwrite(store, call, action) {
        window.done = false;
        window.taskCounter += 1;
        window.collection = [];
        call(action);
        window.starter[taskCounter] = await api.data.getDataValue("segmentation", api.tracing.getCameraPosition());
    }
    async function moveFlycamOverwrite(store, call, action) {
        if (!window.done) {
            call(action);
        }
    }
    api.utils.registerOverwrite("CREATE_NODE", createNodeOverwrite);
    api.utils.registerOverwrite("WK_READY", wKReadyOverwrite);
    api.utils.registerOverwrite("MOVE_FLYCAM", moveFlycamOverwrite);
    api.utils.registerKeyHandler("9", () => {
        api.utils.showToast("info", 'Task finished', 3000);
        api.tracing.finishAndGetNextTask();
    });
    api.utils.registerKeyHandler("7", () => { 
        api.tracing.setCommentForNode("unsolved", api.tracing.getActiveNodeId(), api.tracing.getActiveTreeId());
        api.utils.showToast("warning", 'Task marked as unsolved (mistracing, glia, unsure).', 3000);
        api.tracing.finishAndGetNextTask();
    });
    window.starter[taskCounter] = await api.data.getDataValue("segmentation", api.tracing.getCameraPosition());
})
