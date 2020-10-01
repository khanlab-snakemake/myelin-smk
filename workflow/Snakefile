from os.path import join, split
import pandas as pd

configfile: 'config.yaml'

# Load participants.tsv file
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 

wildcard_constraints:
    subject = "[-a-zA-Z0-9]+"

rule all:
    input: expand('output/myelin_surface/{subject}/lh.myelin.mgh', subject=subjects)


def collect_input(wildcards):
    bids = config['bids']
    subject = '{wildcards.subject}'.format(wildcards=wildcards)

    return { 'T1w': join(bids,'sub-{s}/anat/sub-{s}{f}'.format(s=subject,f=config['scans']['t1w'])),
             'T2w': join(bids,'sub-{s}/anat/sub-{s}{f}'.format(s=subject,f=config['scans']['t2w'])),
             'freesurfer': join(config['fs_deriv'],'{}'.format(subject)) }

# Generate brain mask for more robust boundary-based registration (BBR)

rule get_brain_mask:
    input: unpack(collect_input)
    output: 
        oT2w_brain = 'output/brainmask/{subject}/brain.nii.gz',
        oT2w_brain_mask = 'output/brainmask/{subject}/brain_mask.nii.gz'
    group: 'general'
    singularity: config['neuroglia']        
    shell:
        "out=`basename {output.oT2w_brain} .nii.gz` && "
        "bet {input.T2w} {output.oT2w_brain} -m -R -f 0.25 -v"

# Run FreeSurfer's BBR binary

rule coreg_T2w:
    input: 
        unpack(collect_input),
        T2w_brain = rules.get_brain_mask.output.oT2w_brain
    output:
        oT2w_bbr = 'output/coreg_T2w/{subject}/T2w_reg.dat',
        oT2w = 'output/coreg_T2w/{subject}/T2w_coreg.nii.gz',
        oT1w = 'output/coreg_T2w/{subject}/T1w_coreg.nii.gz',
        obrain = 'output/coreg_T2w/{subject}/brain_mask.nii.gz',        
    params:
        license = config['freesurfer_license'],
        subjects_dir = config['fs_deriv'],
        subject = '{subject}'
    group: 'general'
    threads: 8
    resources:
        mem_mb = 32000        
    singularity: config['freesurfer']
    shell:
        "export FS_LICENSE={params.license} && SUBJECTS_DIR={params.subjects_dir} && "
        "bbregister --s {params.subject} --mov {input.T2w_brain} --reg {output.oT2w_bbr} --init-header --t2 --o {output.oT2w} && "
        "mri_convert {input.T1w} {output.oT1w} -rl $SUBJECTS_DIR/{params.subject}/mri/orig.mgz -nc && "
        "mri_binarize --i $SUBJECTS_DIR/{params.subject}/mri/aparc+aseg.mgz --min 1 --o {output.obrain}"

## Following rules are based on HCP's minimal preprocessing pipeline scripts        
# "PreFreeSurfer/scripts/BiasFieldCorrection_sqrtT1wXT2w.sh"

rule calculate_bias_field:
    input:
        T1w = rules.coreg_T2w.output.oT1w,
        T2w = rules.coreg_T2w.output.oT2w,
        brain = rules.coreg_T2w.output.obrain
    output: 
        obias = 'output/bias_field/{subject}/BiasField.nii.gz',
        oT1w = 'output/bias_field/{subject}/T1w_bc.nii.gz',
        oT2w = 'output/bias_field/{subject}/T2w_bc.nii.gz'
    params:
        sigma = 5
    group: 'general'        
    singularity: config['neuroglia']          
    shell:
        "bash scripts/biascorrection.sh {input.T1w} {input.T2w} {input.brain} {output.obias} {output.oT1w} {output.oT2w} {params.sigma}"
    
    
rule calculate_myelin_map:
    input:
        T1w = rules.calculate_bias_field.output.oT1w,
        T2w = rules.calculate_bias_field.output.oT2w
    output: 'output/myelin_volume/{subject}/T1wDividedByT2w.nii.gz'
    group: 'general'    
    singularity: config['connectome_workbench']  
    shell:
        "wb_command -volume-math 'clamp((T1w / T2w), 0, 100)' {output} -var T1w {input.T1w} -var T2w {input.T2w} -fixnan 0 && "
        "wb_command -volume-palette {output} MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false"


rule map_myelin_to_surface:
    input: rules.calculate_myelin_map.output
    output:
        native = expand('output/myelin_surface/{{subject}}/{hemi}.myelin.{ext}', hemi=['lh','rh'], ext=['mgh']), #,'shape.gii' 
        fsaverage = expand('output/myelin_surface/{{subject}}/{hemi}.myelin.fsaverage.{ext}', hemi=['lh','rh'], ext=['mgh','shape.gii']) 
    params:
        license = config['freesurfer_license'],
        subjects_dir = config['fs_deriv'],
    group: 'general'    
    singularity: config['freesurfer']               
    shell:
        "export FS_LICENSE={params.license} && SUBJECTS_DIR={params.subjects_dir} && "
        "for hemi in lh rh ; do "
            "mri_vol2surf --mov {input} --regheader {wildcards.subject} --cortex "
            "--hemi $hemi --projfrac 0.5 --o `dirname {output.native[0]}`/$hemi.myelin.mgh ; "
            "mri_surf2surf --srcsubject {wildcards.subject} --sval `dirname {output.native[0]}`/$hemi.myelin.mgh "
            "--trgsubject fsaverage --tval `dirname {output.native[0]}`/$hemi.myelin.fsaverage.mgh --hemi $hemi ; "
            ""
            "mri_convert `dirname {output.native[0]}`/$hemi.myelin.mgh `dirname {output.native[0]}`/$hemi.myelin.shape.gii ; "            
            "mri_convert `dirname {output.native[0]}`/$hemi.myelin.fsaverage.mgh `dirname {output.native[0]}`/$hemi.myelin.fsaverage.shape.gii ;"            
        "done"
    
    
#rule surface_mapping:
#    input: unpack(collect_input)
#    output: expand('output/surface_mapping/{{subject}}/{h}.{{subject}}_HCP-MMP.fsaverage.shape.gii', h=hemis, m=modalities)
#    params:
#        subject = '{subject}',
#        subjects_dir = config['subjects_dir'],
#        out_dir = directory('output/surface_mapping/{subject}')
#    run:
#        "bash scripts/fsmapping.sh {params.subject} {input.uni_orig} {input.uni_corr} {input.t1_orig} {input.t1_corr} "
#        "{params.out_dir} {params.subjects_dir}"
