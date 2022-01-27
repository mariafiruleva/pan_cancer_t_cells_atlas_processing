rule integration:
    output: h5="out/{cond}/{cond}/data.h5"
    log: "logs/seurat/{cond}.log"
    params: data=lambda wildcards, output: f'GSE156728/{ind_names[wildcards.cond]}',
            out_dir=lambda wildcards, output: f'out/{wildcards.cond}',
            sample_id=lambda wildcards, output: wildcards.cond,
            annot = config['annot'],
            rscript = config['rscript']
    threads: 4
    shell: """
           {params.rscript} workflow/scripts/integration.R \
           --data {params.data} \
           --out_dir {params.out_dir} \
           --annot {params.annot} \
           --sample_id {params.sample_id} 2> {log}
           """
           
rule integration_by_cd_type:
    output: h5="out/{cd_type}/{cd_type}/data.h5"
    log: "logs/seurat/{cd_type}.log"
    params: data=lambda wildcards, output: ind_by_cd_names[wildcards.cd_type],
            out_dir=lambda wildcards, output: f'out/{wildcards.cd_type}',
            sample_id=lambda wildcards, output: f'{wildcards.cd_type}',
            annot = config['annot'],
            rscript = config['rscript']
    threads: 4
    shell: """
           {params.rscript} workflow/scripts/integrate_multiple_mats.R \
           --data {params.data} \
           --out_dir {params.out_dir} \
           --annot {params.annot} \
           --sample_id {params.sample_id} 2> {log}
           """