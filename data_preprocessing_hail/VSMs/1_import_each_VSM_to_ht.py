import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 

# # ESM1b
# esm1b_ht = hl.import_table(
#         'gs://missense-scoring/esm1b-650m-brandes/proc/*.txt.bgz',
#         types = {
#             'aa_pos':hl.tint,
#             'aa_ref':hl.tstr,
#             'aa_alt':hl.tstr,
#             'esm_score':hl.tfloat,
#             'uniprot':hl.tstr,
#         },
#         force=True,
#         missing='',
#     )
# esm1b_ht = esm1b_ht.rename({'uniprot': 'uniprot_isoform'})
# esm1b_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/esm1b_scores.ht')

# # MisFit
# misfit_ht = hl.import_table(
#         'gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/*.txt.gz',
#         source_file_field='filename',
#         types = {
#             'Uniprot_position':hl.tint,
#             'AA_alt':hl.tstr,
#             'MisFit_D':hl.tfloat,
#             'MisFit_S':hl.tfloat,
#         },
#         force=True,
#         missing='',
#     )
# misfit_ht = misfit_ht.annotate(
#     uniprot_id = misfit_ht.filename.split('/')[-1].split('\.')[0],
# )
# misfit_ht = misfit_ht.rename({'Uniprot_position': 'aa_pos', 'AA_alt': 'aa_alt'})
# misfit_ht = misfit_ht.drop('filename')

# mf_mapping = hl.import_table('gs://missense-scoring/misfit/raw/geneset_s_gene.txt')
# mf_mapping = mf_mapping.key_by('UniprotID')
# misfit_ht = misfit_ht.annotate(
#     enst = mf_mapping[misfit_ht.uniprot_id].TranscriptID, 
#     ensg = mf_mapping[misfit_ht.uniprot_id].GeneID, 
#     ensp = mf_mapping[misfit_ht.uniprot_id].ProteinID, 
#     gene_symbol = mf_mapping[misfit_ht.uniprot_id].Symbol
# )
# misfit_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/misfit_scores.ht')

# # PopEVE
# popeve_raw_ht = hl.import_table(
#     'gs://missense-scoring/popEVE_ukbb_20250312/*.csv', # 2025 popeve scores
#     delimiter=',',
#     source_file_field='filename',
#     types = {
#         'popEVE':hl.tfloat,
#         'popped EVE':hl.tfloat,
#         'popped ESM-1v':hl.tfloat,
#         'EVE':hl.tfloat,
#         'ESM-1v':hl.tfloat,
        
#     },
#     force=True,
#     missing='',
# )

# popeve_ht = popeve_raw_ht.select(
#     aa_pos = hl.int(popeve_raw_ht.mutant[1:-1]),
#     aa_ref = popeve_raw_ht.mutant[0],
#     aa_alt = popeve_raw_ht.mutant[-1],
#     popEVE = popeve_raw_ht['popEVE'], 
#     popped_EVE = popeve_raw_ht['popped EVE'], 
#     popped_ESM_1v = popeve_raw_ht['popped ESM-1v'], 
#     EVE = popeve_raw_ht['EVE'], 
#     ESM_1v = popeve_raw_ht['ESM-1v'], 
#     NP = popeve_raw_ht.filename.split('/')[-1][:-4]
#     )
# popeve_ht.write('gs://missense-scoring/popeve_2025.ht')
# popeve_ht.write('gs://genetics-gym-not-public/Trisha/updated_data_tables/popeve_2025.ht')

# # MPC
# # original: 
# ht = hl.read_table('gs://asc-v17/mpc_gnomad_2.1.1/mpc_grch38_deduped_with_outliers_2024-04-30.ht')
# ht = ht.key_by('locus','alleles')
# ht.write('gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht')

# mpc_ht = hl.read_table('gs://missense-scoring/mpc_grch38_deduped_with_outliers_2024-04-30.ht')

# # RASP
# rasp_table = 'gs://nnfc-fdp-konrad-public/RaSP/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.ht'
# ht = hl.read_table(rasp_table)
# ht = ht.naive_coalesce(50)
# ht = ht.annotate(
#     uniprot_id = ht.uniprot,
#     uniprot_aa_pos = hl.int(ht.variant[1:-1]),
#     uniprot_aa_ref = ht.variant[0],
#     uniprot_aa_alt = ht.variant[-1],
#     rasp_score = ht.score_ml,
# )
# ht = ht.key_by(
#     'uniprot_id',
#     'uniprot_aa_pos',
#     'uniprot_aa_ref',
#     'uniprot_aa_alt',
# )
# ht = ht.select(ht.rasp_score)
# ht.write('gs://missense-scoring/mutation/rasp_scores.ht')

# # AM
# # all isoforms table
# am_table = hl.import_table('gs://dm_alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz', 
#                            delimiter='\t', force_bgz=True, no_header=False, comment='#')
# am_table = am_table.annotate(
#    enst = am_table.transcript_id.split('\.')[0], 
#    aa_pos = hl.int(am_table.protein_variant[1:-1]),
#    aa_ref = am_table.protein_variant[0],
#    aa_alt = am_table.protein_variant[-1],
#    AM = hl.float(am_table.am_pathogenicity), 
#    enst_orig = am_table.transcript_id
# )
# am_table = am_table.select('enst', 'aa_pos', 'aa_ref', 'aa_alt', 'AM', 'enst_orig')
# am_table.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM_isos.ht')

# # # canonical table
# am_table = hl.import_table('gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz', 
#                            delimiter='\t', force_bgz=True, no_header=True, comment='#')
# am_table_fixed = am_table.annotate(
#     locus = hl.locus(am_table.f0, hl.int(am_table.f1), reference_genome='GRCh38'),
#     alleles = [am_table.f2, am_table.f3], 
#     AM_score = hl.float(am_table.f8), 
#     enst_orig = am_table.f6, 
#     consequence = am_table.f9, 
#     genome = am_table.f4
#     )
# am_table = am_table_fixed.annotate(enst = am_table_fixed.enst_orig.split('\.')[0])
# am_table_fixed = am_table_fixed.drop('f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9')
# am_table_fixed_k = am_table_fixed.key_by('locus', 'alleles')
# am_table_fixed_k.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/AM.ht')

# # Revel
# revel_table = hl.import_table('gs://missense-scoring/revel_with_transcript_ids', delimiter=',')
# revel_table = revel_table.filter(~(revel_table.grch38_pos == '.'))
# revel_table = revel_table.annotate(
#     locus = hl.locus('chr' + revel_table.chr, hl.int(revel_table.grch38_pos), reference_genome='GRCh38'),
#     alleles = [revel_table.ref, revel_table.alt], 
#     aa_ref = revel_table.aaref, 
#     aa_alt = revel_table.aaalt, 
#     enst = revel_table.Ensembl_transcriptid.split(';'), 
#     revel = hl.float(revel_table.REVEL)
# )
# revel_table = revel_table.explode('enst')
# revel_table_fixed = revel_table.select('locus', 'alleles', 'enst', 'revel' )
# revel_table_fixed.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/revel_enst.ht')

# ## CPT
# cpt_ht = hl.import_table('gs://missense-scoring/cpt_all/*.csv.gz', 
#                         delimiter=',', force=True, source_file_field='filename')
# cpt_ht = cpt_ht.annotate(
#     aa_pos = hl.int(cpt_ht.mutant[1:-1]),
#     aa_ref = cpt_ht.mutant[0],
#     aa_alt = cpt_ht.mutant[-1],
#     cpt_score = hl.float(cpt_ht.CPT1_score), 
#     uniprot_entry = cpt_ht.filename.split('/')[-1].split('\.')[0]
# )
# cpt_ht = cpt_ht.key_by('aa_pos', 'aa_ref', 'aa_alt')
# cpt_ht.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/cpt_all.ht')

# ## ProteinMPNN
# mpnn_ht = hl.read_table('gs://missense-scoring/mpnn-outputs/AF_total_variants.ht')
# mpnn_ht = mpnn_ht.key_by('uniprot')
# # get isoform mapping from AF name
# af_ht = hl.read_table('gs://genetics-gym-not-public/Emily/af2db_uniprot_final.ht')
# mpnn_ht = mpnn_ht.join(af_ht, how='left')
# mpnn_ht = mpnn_ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')
# mpnn_ht = mpnn_ht.select()
# mpnn_ht.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/protein_mpnn.ht')

# ## PAI3D
# ht = hl.import_csv('gs://missense-scoring/primate_ai3d/PrimateAI-3D_scores.csv.bgz', min_partitions=100)
# ht = ht.key_by(locus=hl.locus(ht.chr, hl.int(ht.pos), 'GRCh38'),
#                alleles=[ht.non_flipped_ref, ht.non_flipped_alt],
#                enst = ht.gene_name.split('\.')[0])
# ht = ht.annotate(score_PAI3D=hl.float(ht.score_PAI3D))
# ht = ht.select(ht.score_PAI3D)
# ht.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/PAI3D.ht')

# ## Polyphen2
# ht = hl.read_table('gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht/')
# ht = ht.explode(ht.vep.transcript_consequences)
# ht = ht.filter(ht.vep.transcript_consequences.transcript_id.startswith('ENST') )
# ht = ht.annotate(
#     polyphen_score = ht.vep.transcript_consequences.polyphen_score, 
#     enst = ht.vep.transcript_consequences.transcript_id, 
#     ensg = ht.vep.transcript_consequences.gene_id, 
#     gene_symbol = ht.ranscript_consequences.gene_symbol)
# ht = ht.select(ht.enst, ht.polyphen_score, ht.ensg, ht.gene_symbol)
# ht = ht.filter(hl.is_defined(ht.polyphen_score))
# ht.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/polyphen.ht')

## CADD
ht = hl.import_table('gs://genetics-gym-not-public/Trisha/whole_genome_SNVs.tsv.gz', 
            filter = '##', delimiter = '\t', force_bgz =True)
ht = ht.select(
    chrom = ht['#Chrom'], 
    pos = hl.int(ht['Pos']),
    ref = ht['Ref'], 
    alt = ht['Alt'], 
    raw_score = hl.float(ht['RawScore']), 
    
)
ht = ht.select(
    locus = hl.locus('chr' + ht['chrom'],ht['pos'], reference_genome='GRCh38'),
    alleles = [ht['ref'], ht['alt']],
    cadd_score = ht['raw_score']
)
ht = ht.key_by('locus', 'alleles')
ht = ht.filter(hl.is_defined(ht.cadd_score))
ht.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/cadd.ht', overwrite=True)


## GPN-MSA
ht = hl.import_table('gs://missense-scoring/GPN-MSA/scores.tsv.bgz', delimiter='\t', force_bgz=True, no_header=True)
ht = ht.annotate(
    chr = 'chr'+ht.f0, 
    pos = hl.int(ht.f1), 
    alleles=[ht.f2, ht.f3], 
    gpn_msa_score = hl.float(ht.f4)
    )
ht = ht.annotate(
    locus = hl.locus(ht.chr, ht.pos, reference_genome='GRCh38'),
)
ht = ht.select(ht.locus, ht.alleles, ht.gpn_msa_score)
ht = ht.key_by('locus', 'alleles')
ht.write('gs://genetics-gym-not-public/Trisha/hail_VSM_tables_updated/gpn_msa.ht')





