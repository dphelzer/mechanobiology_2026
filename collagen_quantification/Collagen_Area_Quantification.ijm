// Collagen I IF area quantification. Fixed intensity-processing params for
// all images; manual steps are spatial only (tissue ROI, artifacts, manual
// patches, post-mask exclusion). Saves final mask per image as *_mask.tif.

var rollingBallRadius = 50;
var medianRadius = 1;
var contrastSaturated = 0.35;

var localThreshMethod = "Phansalkar";
var phansalkarK = 0.5;

var localThreshRadiusFine = 15;
var localThreshRadiiCoarse = newArray(20, 30, 50);

var openIterations = 1;
var closeIterations = 3;

var useHysteresisGrowth = true;
var hysteresisMeanRadius = 150;
var hysteresisStdevRadius = 50;
var hysteresisZThreshold = 0.35;
var hysteresisGrowIterations = 20;

var usePerImagePatchCalibration = true;
var patchThreshMarginStdevs = 2;
var manualPatchThreshValue = 22; // fallback if usePerImagePatchCalibration=false

var useIlluminationFlatten = true;
var flattenSigmaUm = 2500; // must be >> size of brightness-shifted region

var minParticleSizeUm2 = 5;

var debugPauseAtMask = true; // also enables post-mask exclusion loop

var flattenSigmaPx = 0;
var manualPatchRoiIndex = -1;
var activePatchThreshValue = manualPatchThreshValue;

var outputDir = "";
var qcMean = 0;
var qcStd = 0;

// ============================================================================
macro "Collagen Area Quantification" {

    if (useIlluminationFlatten) {
        flattenStatus = "ENABLED, sigma=" + flattenSigmaUm + " um";
    } else {
        flattenStatus = "DISABLED";
    }
    if (usePerImagePatchCalibration) {
        patchStatus = "PER-IMAGE (neg. mean + " + patchThreshMarginStdevs + "*stdev)";
    } else {
        patchStatus = "FIXED, cutoff=" + manualPatchThreshValue;
    }
    if (useHysteresisGrowth) {
        hysStatus = "ENABLED (z-score, mean r=" + hysteresisMeanRadius + ", stdev r=" + hysteresisStdevRadius + ", threshold=" + hysteresisZThreshold + ", max " + hysteresisGrowIterations + " steps)";
    } else {
        hysStatus = "DISABLED";
    }

    Dialog.create("Locked parameters check");
    Dialog.addMessage(
        "These parameters are fixed for ALL images in this run:\n\n" +
        "Rolling ball radius: " + rollingBallRadius + " px\n" +
        "Median despeckle radius: " + medianRadius + " px\n" +
        "Contrast normalization: saturated=" + contrastSaturated + "%\n" +
        "Local threshold, fine: " + localThreshMethod + " radius " + localThreshRadiusFine + ", k=" + phansalkarK + "\n" +
        "Local threshold, coarse radii: " + String.join(localThreshRadiiCoarse, ",") + "\n" +
        "Manually-flagged patch threshold: " + patchStatus + "\n" +
        "Illumination flattening: " + flattenStatus + "\n" +
        "Open iterations: " + openIterations + " / Close iterations: " + closeIterations + "\n" +
        "Hysteresis growth: " + hysStatus + "\n" +
        "Min particle size: " + minParticleSizeUm2 + " um^2\n\n" +
        "Do NOT change these between samples in the same experiment.\n" +
        "Click OK to confirm and continue.");
    Dialog.show();

    // Must be set before any thresholding, or Convert to Mask silently inverts.
    run("Options...", "iterations=1 count=1 black");

    run("Set Measurements...", "area area_fraction limit redirect=None decimal=3");
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    if (isOpen("Summary")) { selectWindow("Summary"); run("Close"); }

    outputDir = getDirectory("Choose folder to save masks, ROIs, and results CSV");

    continueLoop = true;
    while (continueLoop) {
        processOneImage();
        continueLoop = getBoolean("Process another image?");
    }

    if (nResults > 0) {
        saveAs("Results", outputDir + "Collagen_Area_Results.csv");
        print("Saved results to " + outputDir + "Collagen_Area_Results.csv");
    }
    print("Done.");
}

// ============================================================================
function processOneImage() {

    waitForUser("SELECT IMAGE",
        "Click the raw image window to make it active.\nClick OK to continue.");
    origTitle = getTitle();

    getVoxelSize(vw, vh, vd, unit);
    if (unit == "pixel" || unit == "" || unit == "pixels") {
        showMessage("No spatial calibration detected",
            "This image has no real-world scale set. Areas will NOT be\n" +
            "comparable across samples unless every image is calibrated in\n" +
            "physical units matching acquisition metadata.\n\n" +
            "The Properties dialog will open now -- enter the pixel size.");
        run("Properties...");
        getVoxelSize(vw, vh, vd, unit);
    }

    run("Duplicate...", "title=[" + origTitle + "_work]");
    workTitle = getTitle();

    analysisRoiIndex = defineAnalysisRoi();
    computeFlattenSigmaPx();

    grayTitle = preprocessToNormalized8bit(analysisRoiIndex);

    // Hybrid mask built on the full unclipped image; clip happens later.
    selectWindow(grayTitle);
    run("Duplicate...", "title=local_mask_fine");
    run("Auto Local Threshold", "method=" + localThreshMethod +
        " radius=" + localThreshRadiusFine + " parameter_1=" + phansalkarK + " parameter_2=0 white");
    rename("accum_mask");

    for (r = 0; r < localThreshRadiiCoarse.length; r++) {
        radius = localThreshRadiiCoarse[r];
        selectWindow(grayTitle);
        run("Duplicate...", "title=temp_mask");
        run("Auto Local Threshold", "method=" + localThreshMethod +
            " radius=" + radius + " parameter_1=" + phansalkarK + " parameter_2=0 white");
        imageCalculator("OR create", "accum_mask", "temp_mask");
        close("accum_mask");
        close("temp_mask");
        rename("accum_mask");
    }

    // Hysteresis growth, seeded from local threshold only; patch mask added after.
    if (useHysteresisGrowth) {
        selectWindow(grayTitle);
        run("Duplicate...", "title=local_mean_img");
        run("Mean...", "radius=" + hysteresisMeanRadius);

        selectWindow(grayTitle);
        run("Duplicate...", "title=local_var_img");
        run("32-bit");
        run("Variance...", "radius=" + hysteresisStdevRadius);
        run("Square Root");
        run("Add...", "value=1");

        imageCalculator("Subtract create 32-bit", grayTitle, "local_mean_img");
        close("local_mean_img");
        rename("elevation_img");

        imageCalculator("Divide create 32-bit", "elevation_img", "local_var_img");
        close("elevation_img");
        close("local_var_img");
        rename("loose_mask");
        setThreshold(hysteresisZThreshold, 1e30);
        run("Options...", "iterations=1 count=1 black");
        run("Convert to Mask");

        // Clip seed to tissue-minus-artifacts ROI, and re-clip every growth
        // iteration since Dilate ignores ROI boundaries.
        selectWindow("accum_mask");
        run("Duplicate...", "title=grown_mask");
        roiManager("Select", analysisRoiIndex);
        run("Make Inverse");
        setColor(0);
        fill();
        run("Select None");
        run("Options...", "iterations=1 count=1 black");
        for (g = 0; g < hysteresisGrowIterations; g++) {
            selectWindow("grown_mask");
            run("Dilate");
            imageCalculator("AND", "grown_mask", "loose_mask");
            roiManager("Select", analysisRoiIndex);
            run("Make Inverse");
            setColor(0);
            fill();
            run("Select None");
        }
        // Replace accum_mask with grown_mask -- do not union with original.
        close("accum_mask");
        selectWindow("grown_mask");
        rename("accum_mask");
        close("loose_mask");
    }

    // Manually-flagged patch mask, added after hysteresis growth.
    if (manualPatchRoiIndex >= 0) {
        if (usePerImagePatchCalibration) {
            patchMaskAccepted = false;
            while (!patchMaskAccepted) {
                selectWindow(grayTitle);
                run("Select None");
                negStartIdx = roiManager("count");
                waitForUser("SAMPLE NEGATIVE TISSUE (this image)",
                    "Draw a selection entirely INSIDE a region you are confident\n" +
                    "has NO collagen, in THIS image, then press T to add it to\n" +
                    "the ROI Manager. Repeat for 2-4 different spots (recommended\n" +
                    "-- pooling averages out regional background drift).\n\n" +
                    "When you've added all of them, click OK. You can zoom/pan\n" +
                    "freely throughout this step.");
                nNegSamples = roiManager("count") - negStartIdx;
                if (nNegSamples == 0) {
                    exit("No negative regions were added. Rerun and press T after drawing each one.");
                }
                negIndices = newArray(nNegSamples);
                for (i = 0; i < nNegSamples; i++) negIndices[i] = negStartIdx + i;
                roiManager("Select", negIndices);
                roiManager("OR");
                getStatistics(calNegArea, calNegMean, calNegMin, calNegMax, calNegStd);
                run("Select None");
                activePatchThreshValue = calNegMean + patchThreshMarginStdevs * calNegStd;
                print(origTitle + ": per-image patch threshold = " +
                      d2s(activePatchThreshValue, 1) + " (pooled from " +
                      nNegSamples + " samples, mean " + d2s(calNegMean, 1) +
                      " + " + patchThreshMarginStdevs + "*stdev " + d2s(calNegStd, 1) + ")");

                generatePatchMask(activePatchThreshValue);

                selectWindow("manual_patch_mask");
                patchMaskAccepted = getBoolean(
                    "Patch mask (threshold=" + d2s(activePatchThreshValue, 1) +
                    ") is on screen. Does it look correct?\n\n" +
                    "No = discard and resample the negative control regions.");
                if (!patchMaskAccepted) {
                    close("manual_patch_mask");
                }
            }
        } else {
            activePatchThreshValue = manualPatchThreshValue;
            generatePatchMask(activePatchThreshValue);
        }
        imageCalculator("OR", "accum_mask", "manual_patch_mask");
        close("manual_patch_mask");
    }

    close(grayTitle);
    selectWindow("accum_mask");
    rename("combined_mask");

    if (debugPauseAtMask) {
        waitForUser("INSPECT RAW MASK (pre-correction)",
            "Raw combined mask, BEFORE clipping/Open/Close. Edge brightness\n" +
            "and isolated spots may still be visible here -- judge whether\n" +
            "corrections worked at the FINAL MASK checkpoint later.\n\n" +
            "Click OK to continue.");
    }

    roiManager("Select", analysisRoiIndex);
    run("Make Inverse");
    setColor(0);
    fill();
    run("Select None");

    run("Options...", "iterations=" + openIterations + " count=1 black");
    run("Open");

    run("Options...", "iterations=" + closeIterations + " count=1 black");
    run("Close-");

    selectWindow("combined_mask");
    roiManager("Select", analysisRoiIndex);
    run("Analyze Particles...", "size=" + minParticleSizeUm2 + "-Infinity show=Masks summarize");
    finalMaskTitle = getTitle(); // this is what gets measured and saved

    finalRoiIndex = analysisRoiIndex;

    if (debugPauseAtMask) {
        satisfiedWithMask = false;
        while (!satisfiedWithMask) {
            wantsExclusion = getBoolean(
                "Final mask (after clipping, Open, Close, size filtering) is\n" +
                "on screen -- this is exactly what gets measured and saved.\n\n" +
                "Looking at it, are there specific regions you want to exclude\n" +
                "(e.g. residual bright-edge false positives)? Click Yes to mark\n" +
                "them now, or No to accept this mask as-is.");

            if (!wantsExclusion) {
                satisfiedWithMask = true;
            } else {
                postHocStartIdx = roiManager("count");
                waitForUser("MARK REGIONS TO EXCLUDE",
                    "Draw a selection around a region to exclude, press T to add\n" +
                    "it to the ROI Manager. Repeat for each additional region.\n\n" +
                    "When you've added all of them, click OK. You can zoom/pan\n" +
                    "freely throughout this step.");
                nPostHoc = roiManager("count") - postHocStartIdx;
                if (nPostHoc == 0) {
                    exit("No exclusion regions were added. Rerun and press T after drawing each one.");
                }
                postHocIndices = newArray(nPostHoc);
                for (i = 0; i < nPostHoc; i++) postHocIndices[i] = postHocStartIdx + i;
                roiManager("Select", postHocIndices);
                roiManager("OR");
                roiManager("Add");
                postHocUnionIndex = roiManager("count") - 1;

                // Clip exclusion to current analysis region, then subtract
                roiManager("Select", newArray(finalRoiIndex, postHocUnionIndex));
                roiManager("AND");
                if (selectionType() != -1) {
                    roiManager("Add");
                    postHocClipIndex = roiManager("count") - 1;
                    roiManager("Select", newArray(finalRoiIndex, postHocClipIndex));
                    roiManager("XOR");
                    roiManager("Add");
                    finalRoiIndex = roiManager("count") - 1;
                }

                selectWindow("combined_mask");
                roiManager("Select", postHocUnionIndex);
                setColor(0);
                fill();
                run("Select None");

                close(finalMaskTitle);
                selectWindow("combined_mask");
                roiManager("Select", finalRoiIndex);
                run("Analyze Particles...", "size=" + minParticleSizeUm2 + "-Infinity show=Masks summarize");
                finalMaskTitle = getTitle();
            }
        }
    }

    close("combined_mask");

    // --- Export tissue (analysis region) mask for downstream analysis ---
    selectWindow(finalMaskTitle);
    w = getWidth();
    h = getHeight();
    newImage("tissue_mask", "8-bit black", w, h, 1);
    setVoxelSize(vw, vh, vd, unit);
    roiManager("Select", finalRoiIndex);
    setColor(255);
    fill();
    // Clip to the original tissue outline (see note below)
    roiManager("Select", 0);
    run("Make Inverse");
    setColor(0);
    fill();
    run("Select None");
    tissueSafeName = replace(origTitle, "[^A-Za-z0-9_\\-\\.]", "_");
    saveAs("Tiff", outputDir + tissueSafeName + "_tissue.tif");
    close();

    // --------------------------------------------------------------------

    sumRows     = Table.size("Summary");
    totalArea   = Table.get("Total Area", sumRows - 1, "Summary");
    percentArea = Table.get("%Area", sumRows - 1, "Summary");
    count       = Table.get("Count", sumRows - 1, "Summary");

    row = nResults;
    setResult("Image", row, origTitle);
    setResult("Collagen Area (calibrated)", row, totalArea);
    setResult("Percent Area (%)", row, percentArea);
    setResult("Particle Count", row, count);
    setResult("Pixel width", row, vw);
    setResult("Unit", row, unit);
    setResult("QC PreNorm Mean Intensity", row, qcMean);
    setResult("QC PreNorm StdDev", row, qcStd);
    updateResults();

    selectWindow(finalMaskTitle);
    safeName = replace(origTitle, "[^A-Za-z0-9_\\-\\.]", "_");
    saveAs("Tiff", outputDir + safeName + "_mask.tif");
    roiManager("Deselect");
    roiManager("Save", outputDir + safeName + "_ROIs.zip");

    print(origTitle + ":  Area=" + totalArea + "  %Area=" + percentArea +
          "  Count=" + count);

    close();
}

// Thresholds grayTitle at threshValue, restricted to manualPatchRoiIndex,
// leaving "manual_patch_mask" as the active window.
function generatePatchMask(threshValue) {
    selectWindow(grayTitle);
    run("Duplicate...", "title=manual_patch_mask");
    setThreshold(threshValue, 255);
    run("Options...", "iterations=1 count=1 black");
    run("Convert to Mask");
    roiManager("Select", manualPatchRoiIndex);
    run("Make Inverse");
    setColor(0);
    fill();
    run("Select None");
}

// ============================================================================
// Draws tissue ROI, then optionally marks artifact (excluded) and diffuse-
// patch (kept, thresholded) regions. Returns ROI Manager index of the
// tissue-minus-artifacts composite. Sets global manualPatchRoiIndex (-1 if
// none). All intermediate ROIs stay in the Manager for the audit trail.
// ============================================================================
function defineAnalysisRoi() {
    roiManager("Reset");
    waitForUser("DEFINE TISSUE ROI",
        "Draw a polygon ROI around the tissue region, press T to add it,\n" +
        "then click OK.\n\n" +
        "If a locally bright region (e.g. the typical IF edge effect) needs\n" +
        "excluding, you can inset this outline, mark it as an artifact in the\n" +
        "next step, or exclude it later after seeing the mask -- whichever is\n" +
        "easiest for that specific case.");
    if (roiManager("count") == 0) {
        exit("No tissue ROI was added to the ROI Manager.");
    }

    artifactStartIdx = roiManager("count");
    addMore = getBoolean(
        "Are there artifacts (folds, bubbles, dust) inside the tissue that\n" +
        "need excluding from analysis?");
    if (addMore) {
        waitForUser("DEFINE ARTIFACT REGIONS",
            "Draw a selection around an artifact region, press T to add it to\n" +
            "the ROI Manager. Repeat for each additional artifact. This does\n" +
            "NOT alter the image -- it only marks regions for exclusion.\n\n" +
            "When you've added all of them, click OK. You can zoom/pan freely\n" +
            "throughout this step.");
    }
    nArtifacts = roiManager("count") - artifactStartIdx;

    if (nArtifacts == 0) {
        analysisIndex = 0;
    } else {
        artifactIndices = newArray(nArtifacts);
        for (i = 0; i < nArtifacts; i++) artifactIndices[i] = artifactStartIdx + i;
        roiManager("Select", artifactIndices);
        roiManager("OR");
        roiManager("Add");
        artifactUnionIndex = roiManager("count") - 1;

        // Clip artifacts to tissue so XOR acts as true subtraction
        roiManager("Select", newArray(0, artifactUnionIndex));
        roiManager("AND");
        if (selectionType() == -1) {
            analysisIndex = 0; // artifacts don't overlap tissue: nothing to remove
        } else {
            roiManager("Add");
            artifactClipIndex = roiManager("count") - 1;
            roiManager("Select", newArray(0, artifactClipIndex));
            roiManager("XOR");
            roiManager("Add");
            analysisIndex = roiManager("count") - 1;
        }
    }

    manualPatchRoiIndex = -1;
    addPatch = getBoolean(
        "Are there large diffuse collagen patches you can identify visually\n" +
        "that automated detection can't handle (e.g. touching the tissue\n" +
        "edge)? You'll mark them now; a fixed threshold will be applied\n" +
        "within them.");
    patchStartIdx = roiManager("count");
    if (addPatch) {
        waitForUser("MARK DIFFUSE PATCH REGIONS",
            "Draw a selection around a diffuse patch region, press T to add\n" +
            "it. Repeat for each additional patch region.\n\n" +
            "When you've added all of them, click OK. You can zoom/pan freely\n" +
            "throughout this step.");
    }
    nPatches = roiManager("count") - patchStartIdx;
    if (nPatches > 0) {
        patchIndices = newArray(nPatches);
        for (i = 0; i < nPatches; i++) patchIndices[i] = patchStartIdx + i;
        roiManager("Select", patchIndices);
        roiManager("OR");
        roiManager("Add");
        manualPatchRoiIndex = roiManager("count") - 1;
    }

    run("Select None");
    return analysisIndex;
}

function computeFlattenSigmaPx() {
    flattenSigmaPx = flattenSigmaUm / vw;
}

// Returns the value at loPercentile (0-1) of the current selection's pixel
// distribution -- used to find the dim floor of a diffuse region without
// isolating it by hand.
function getPercentileOfSelection(loPercentile) {
    getStatistics(selArea, selMean, selMin, selMax, selStd);
    nBins = 256;
    getHistogram(values, counts, nBins, selMin, selMax);
    total = 0;
    for (b = 0; b < nBins; b++) total += counts[b];
    target = total * loPercentile;
    running = 0;
    for (b = 0; b < nBins; b++) {
        running += counts[b];
        if (running >= target) {
            return values[b];
        }
    }
    return selMax;
}

// Shared fixed preprocessing: background subtraction -> despeckle ->
// illumination flattening (if enabled) -> QC stats -> contrast normalization
// -> fixed 8-bit conversion. Sets qcMean/qcStd. Returns the result title.
function preprocessToNormalized8bit(roiIndex) {
    run("Subtract Background...", "rolling=" + rollingBallRadius + " sliding");
    run("Median...", "radius=" + medianRadius);

    if (useIlluminationFlatten) {
        run("32-bit");
        rawTitle = getTitle();
        run("Duplicate...", "title=flatten_bg");
        run("Gaussian Blur...", "sigma=" + flattenSigmaPx);
        run("Add...", "value=1");
        imageCalculator("Divide create 32-bit", rawTitle, "flatten_bg");
        close(rawTitle);
        close("flatten_bg");
        run("Multiply...", "value=1000");
        setMinAndMax(0, 65535); // fixed range before 16-bit, avoid auto-scale crushing low values
        run("16-bit");
    }

    roiManager("Select", roiIndex);
    getStatistics(qcArea, qcMean, qcMin, qcMax, qcStd);

    run("Enhance Contrast...", "saturated=" + contrastSaturated + " normalize");
    run("Select None");

    bd = bitDepth();
    if (bd == 8) { theoreticalMax = 255; }
    else if (bd == 16) { theoreticalMax = 65535; }
    else { theoreticalMax = 65535; }
    setMinAndMax(0, theoreticalMax);
    run("8-bit");
    return getTitle();
}

// ============================================================================
// CHECK INTENSITY SEPARATION (diagnostic) -- validates the per-image patch
// threshold formula, fixed fallback, and hysteresis z-threshold.
// ============================================================================
macro "Check Intensity Separation" {
    run("Options...", "iterations=1 count=1 black");

    waitForUser("SELECT IMAGE", "Click a representative image. Click OK.");
    origTitle = getTitle();

    getVoxelSize(vw, vh, vd, unit);
    if (unit == "pixel" || unit == "" || unit == "pixels") {
        showMessage("No spatial calibration detected",
            "Set the image scale via Image > Properties before checking.");
        run("Properties...");
    }

    run("Duplicate...", "title=[" + origTitle + "_calib]");

    analysisRoiIndex = defineAnalysisRoi();
    computeFlattenSigmaPx();
    grayTitle = preprocessToNormalized8bit(analysisRoiIndex);

    run("Select None");
    waitForUser("SAMPLE DIFFUSE COLLAGEN",
        "Draw a selection around the WHOLE diffuse/thick ECM region --\n" +
        "bright and dim parts together. Don't try to isolate just a dim\n" +
        "sub-region -- if it's a smooth gradient with no clean boundary,\n" +
        "sampling the whole patch and using a low percentile of its values\n" +
        "(done automatically below) finds the dim end without needing you\n" +
        "to isolate it. For large patches, sample well away from the tissue\n" +
        "edge. Click OK.");
    getStatistics(diffArea, diffMean, diffMin, diffMax, diffStd);
    roiManager("Add");
    diffRoiIndex = roiManager("count") - 1;

    run("Select None");
    negStartIdx = roiManager("count");
    waitForUser("SAMPLE NEGATIVE TISSUE",
        "Draw a selection entirely INSIDE a region you are confident has NO\n" +
        "collagen, then press T to add it to the ROI Manager. Repeat for\n" +
        "2-4 different spots (recommended, matching the main pipeline).\n\n" +
        "When you've added all of them, click OK. You can zoom/pan freely\n" +
        "throughout this step.");
    nNegSamples = roiManager("count") - negStartIdx;
    if (nNegSamples == 0) {
        exit("No negative regions were added. Rerun and press T after drawing each one.");
    }
    negIndices = newArray(nNegSamples);
    for (i = 0; i < nNegSamples; i++) negIndices[i] = negStartIdx + i;
    roiManager("Select", negIndices);
    roiManager("OR");
    getStatistics(negArea, negMean, negMin, negMax, negStd);
    roiManager("Add");
    negRoiIndex = roiManager("count") - 1;
    run("Select None");

    separation = (diffMean - negMean) / ((diffStd + negStd) / 2);
    negCeiling = negMean + 2 * negStd;
    formulaValue = negMean + patchThreshMarginStdevs * negStd;

    print("=== INTENSITY SEPARATION CHECK: " + origTitle + " ===");
    print("Diffuse region -- mean: " + d2s(diffMean, 1) + "  stdev: " + d2s(diffStd, 1));
    print("Negative region -- mean: " + d2s(negMean, 1) + "  stdev: " + d2s(negStd, 1));
    print("Separation (mean diff / avg stdev): " + d2s(separation, 2));
    print("");
    print("Per-image formula (mean + " + patchThreshMarginStdevs + "*stdev): " + d2s(formulaValue, 1));
    if (formulaValue > negCeiling && formulaValue < diffMean) {
        print("Sits inside this image's gap.");
    } else if (formulaValue >= diffMean) {
        print("WARNING: exceeds the diffuse floor -- may undercount here.");
    } else {
        print("Note: below the negative ceiling -- check if this sample was noisy.");
    }
    print("");
    print("Fallback fixed value (manualPatchThreshValue=" + manualPatchThreshValue + "):");
    if (negCeiling < diffMean) {
        print("Gap: " + d2s(negCeiling, 1) + " to " + d2s(diffMean, 1));
        if (manualPatchThreshValue > negCeiling && manualPatchThreshValue < diffMean) {
            print("Sits inside this image's gap.");
        } else {
            print("WARNING: does NOT sit inside this image's gap.");
        }
    } else {
        print("NO GAP: negative ceiling meets or exceeds the diffuse mean.");
    }
    print("=================================================================");

    selectWindow(grayTitle);
    run("Duplicate...", "title=hyst_local_mean_img");
    run("Mean...", "radius=" + hysteresisMeanRadius);

    selectWindow(grayTitle);
    run("Duplicate...", "title=hyst_local_var_img");
    run("32-bit");
    run("Variance...", "radius=" + hysteresisStdevRadius);
    run("Square Root");
    run("Add...", "value=1");

    imageCalculator("Subtract create 32-bit", grayTitle, "hyst_local_mean_img");
    close("hyst_local_mean_img");
    rename("hyst_elevation_img");
    imageCalculator("Divide create 32-bit", "hyst_elevation_img", "hyst_local_var_img");
    close("hyst_elevation_img");
    close("hyst_local_var_img");
    rename("hyst_z_img");

    roiManager("Select", diffRoiIndex);
    diffZFloor = getPercentileOfSelection(0.10);
    getStatistics(diffZArea, diffZMean, diffZMin, diffZMax, diffZStd);
    roiManager("Select", negRoiIndex);
    getStatistics(negZArea, negZMean, negZMin, negZMax, negZStd);
    run("Select None");
    close("hyst_z_img");

    hystNegCeiling = negZMean + 2 * negZStd;
    hystSeparation = (diffZFloor - negZMean) / ((diffZStd + negZStd) / 2);

    print("=== HYSTERESIS Z-THRESHOLD CALIBRATION: " + origTitle + " ===");
    print("(local z-score, mean radius=" + hysteresisMeanRadius + ", stdev radius=" + hysteresisStdevRadius + ")");
    print("Whole diffuse/ECM sample -- mean: " + d2s(diffZMean, 2) +
          "  stdev: " + d2s(diffZStd, 2) + "  10th percentile (dim floor): " + d2s(diffZFloor, 2));
    print("Negative sample          -- mean: " + d2s(negZMean, 2) + "  stdev: " + d2s(negZStd, 2));
    print("Separation (dim floor vs negative): " + d2s(hystSeparation, 2));
    print("Negative-tissue ceiling (mean + 2*stdev): " + d2s(hystNegCeiling, 2));
    print("Dim floor of diffuse/ECM sample (10th percentile): " + d2s(diffZFloor, 2));
    if (hystNegCeiling < diffZFloor) {
        hystGapMid = (hystNegCeiling + diffZFloor) / 2;
        print("Gap: " + d2s(hystNegCeiling, 2) + " to " + d2s(diffZFloor, 2) +
              " (midpoint " + d2s(hystGapMid, 2) + ").");
        print("Current hysteresisZThreshold (" + hysteresisZThreshold + "):");
        if (hysteresisZThreshold > hystNegCeiling && hysteresisZThreshold < diffZFloor) {
            print("Sits inside this image's gap.");
        } else {
            print("WARNING: does NOT sit inside this image's gap.");
        }
    } else {
        print("NO GAP: negative ceiling meets or exceeds the dim floor -- can be a");
        print("genuine, real limit rather than a calibration error.");
    }
    print("Repeat on several images AND, if possible, several different ECM");
    print("regions -- pick the z-threshold that sits inside every gap.");
    print("=================================================================");

    close();
}
