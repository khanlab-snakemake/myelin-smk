# parse arguments
T1wImage=$1
T2wImage=$2

mask=$3

OutputBiasField=$4
OutputT1wRestoredImage=$5
OutputT2wRestoredImage=$6

BiasFieldSmoothingSigma=$7

# default parameters
Factor="0.5" #Leave this at 0.5 for now it is the number of standard deviations below the mean to threshold the non-brain tissues at
WD=`dirname $OutputBiasField`

T1wImageBrain=$WD/T1w_brain.nii.gz
T2wImageBrain=$WD/T2w_brain.nii.gz

########################################## DO WORK ########################################## 

# Generate brainmask based on T2w image, ie. binarize everything non-zero
fslmaths $T1wImage -mul $mask $T1wImageBrain
fslmaths $T2wImage -mul $mask $T2wImageBrain

# Form sqrt(T1w*T2w), mask this and normalise by the mean
echo " --> Forming sqrt(T1w*T2w), masking this and normalising by the mean"
fslmaths $T1wImage -mul $T2wImage -abs -sqrt $WD/T1wmulT2w.nii.gz -odt float
fslmaths $WD/T1wmulT2w.nii.gz -mas $T1wImageBrain $WD/T1wmulT2w_brain.nii.gz
meanbrainval=`${FSLDIR}/bin/fslstats $WD/T1wmulT2w_brain.nii.gz -M`
fslmaths $WD/T1wmulT2w_brain.nii.gz -div $meanbrainval $WD/T1wmulT2w_brain_norm.nii.gz

# Smooth the normalised sqrt image, using within-mask smoothing : s(Mask*X)/s(Mask)
echo " --> Smoothing the normalised sqrt image, using within-mask smoothing"
fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -bin -s $BiasFieldSmoothingSigma $WD/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz
fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -s $BiasFieldSmoothingSigma -div $WD/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz $WD/T1wmulT2w_brain_norm_s${BiasFieldSmoothingSigma}.nii.gz

# Divide normalised sqrt image by smoothed version (to do simple bias correction)
echo " --> Dividing normalised sqrt image by smoothed version"
fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -div $WD/T1wmulT2w_brain_norm_s$BiasFieldSmoothingSigma.nii.gz $WD/T1wmulT2w_brain_norm_modulate.nii.gz

# Create a mask using a threshold at Mean - 0.5*Stddev, with filling of holes to remove any non-grey/white tissue.
echo " --> Creating a mask and filling holes"
STD=`fslstats $WD/T1wmulT2w_brain_norm_modulate.nii.gz -S`
echo $STD
MEAN=`fslstats $WD/T1wmulT2w_brain_norm_modulate.nii.gz -M`
echo $MEAN
Lower=`echo "$MEAN - ($STD * $Factor)" | bc -l`
echo $Lower
fslmaths $WD/T1wmulT2w_brain_norm_modulate -thr $Lower -bin -ero -mul 255 $WD/T1wmulT2w_brain_norm_modulate_mask
#wb_command -volume-remove-islands $WD/T1wmulT2w_brain_norm_modulate_mask.nii.gz $WD/T1wmulT2w_brain_norm_modulate_mask.nii.gz

# Extrapolate normalised sqrt image from mask region out to whole FOV
echo " --> Extrapolating normalised sqrt image from mask region out to whole FOV"
fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -mas $WD/T1wmulT2w_brain_norm_modulate_mask.nii.gz -dilall $WD/bias_raw.nii.gz -odt float
fslmaths $WD/bias_raw.nii.gz -s $BiasFieldSmoothingSigma $OutputBiasField

# Use bias field output to create corrected images
echo " --> Using bias field output to create corrected images"
#fslmaths $T1wImage -div $OutputBiasField -mas $T1wImageBrain $OutputT1wRestoredBrainImage -odt float
fslmaths $T1wImage -div $OutputBiasField -mas $mask $OutputT1wRestoredImage -odt float
#fslmaths $T2wImage -div $OutputBiasField -mas $T1wImageBrain $OutputT2wRestoredBrainImage -odt float
fslmaths $T2wImage -div $OutputBiasField -mas $mask $OutputT2wRestoredImage -odt float
fslcpgeom $OutputT1wRestoredImage $OutputT2wRestoredImage

echo "---> Finished Bias Field Correction"
