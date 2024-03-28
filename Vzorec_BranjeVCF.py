# from flask import Flask
# from flask_appbuilder import SQLA
# app = Flask(__name__)
# app.config.from_object("config")
# db = SQLA(app)
from app.celery_worker import db

import pandas as pd
from vcf_parser import VCFParser
from sqlalchemy.sql.expression import bindparam
from app.models import Patient, Variants
from app.celery_worker import celery
import logging

# Set the root logger level (optional, defaults to logging.WARNING)
logging.getLogger().setLevel(logging.ERROR)

# patient_identifier="RD401"

@celery.task
def updateVCFdata(patient_identifier):
    current_patient = db.session.query(Patient).filter(Patient.patient_identifier==patient_identifier).first()

    # Define the selected info fields you want to extract
    # selected_info_fields = ["clinvar_classification", "TopMed_AC", "TopMed_AF"]  # Replace with the actual field names you want
    selected_info_fields = ["Regeneron_AC", "Regeneron_AF", "Regeneron_AN", 
                            "TopMed_AC", "TopMed_AF", "TopMed_AN", "TopMed_Het", "TopMed_Hom", 
                            "GnomAD4_exomes_AC", "GnomAD4_exomes_AF", "GnomAD4_exomes_AN", "GnomAD4_exomes_Hom",
                            "GnomAD4_genomes_AC", "GnomAD4_genomes_AF", "GnomAD4_genomes_AN", "GnomAD4_genomes_Hom",
                            "pext_score", "RMC_OE", "RMC_CHISQ",
                            "GERP_RS",
                            "gnomAD_pLi", "gnomAD_misz",
                            "metaDome",
                            "blacklisted_regions"]  # Replace with the actual field names you want

    # Initialize empty lists to store the extracted data
    chromosomes = []
    positions = []
    ref_alleles = []
    alt_alleles = []
    variant_unique_id = []

    info_data = {field: [] for field in selected_info_fields}

    my_parser = VCFParser(infile=current_patient.has_vcf_anno, split_variants=True, check_info=False)
    for variant in my_parser:
        # if "clinvar.CLNSIG" in variant["info_dict"]:
        #     #print(variant["info_dict"]["clinvar.CLNSIG"])

        chromosomes.append(variant["CHROM"])
        positions.append(variant["POS"])
        ref_alleles.append(variant["REF"])
        alt_alleles.append(variant["ALT"])
        variant_unique_id.append("-".join([variant["CHROM"], variant["POS"], variant["REF"], variant["ALT"]]))

        for field in selected_info_fields:
            if field in variant["info_dict"]:
                info_data[field].append(variant["info_dict"][field][0])
            else:
                info_data[field].append(None)

    # Create a DataFrame from the extracted data
    data = {
        'Chromosome': chromosomes,
        'Position': positions,
        'Ref_Allele': ref_alleles,
        'Alt_Allele': alt_alleles,
        'variant_unique_id': variant_unique_id,
    }
    data.update(info_data)
    df = pd.DataFrame(data)

    # Fix pext score column
    df = df.where(pd.notna(df), None)
    df.replace('nan', None, inplace=True)

    # Do not import rows that have all None values - temporarily
    df = df.dropna(subset=selected_info_fields, how='all')

    # Get pathogenicity dict for ClinVar
    # query_result = db.session.query(ClinSig.id, ClinSig.english_name).all()
    # Convert the query result to a list of dictionaries
    # clinvar_map = {row.english_name:row.id for row in query_result}
    # Not working yet
    # df["clinvar.CLNSIG"].map(clinvar_map)
    # df["clinvar_classification"] = df["clinvar.CLNSIG"]

    # Display the DataFrame
    # print(df)

    ########################################################
    # Code for import of variants into the database
    ########################################################
    # Prepare models
    v_model = Variants

    # Prepare bindparams
    fields = ["variant_unique_id"]+selected_info_fields
    v_bindparams = {}
    for value in fields:
        v_bindparams[value]=bindparam(value)

    # Prepare a dict to update variants
    variant_update_dict = df[["variant_unique_id"]+selected_info_fields].to_dict(orient="records")
    print(f'Will update {str(len(variant_update_dict))} variants.') if len(variant_update_dict)>0 else print(f'Update dict empty.')

    # Perform the update
    variant_update_stmt = v_model.__table__.update().\
        where(v_model.__table__.c.variant_unique_id == bindparam('variant_unique_id')).\
        values(v_bindparams)

    batch_size = 1000
    total_dicts = len(variant_update_dict)

    if total_dicts > 0:
        num_batches = (total_dicts + batch_size - 1) // batch_size

        for i in range(num_batches):
            start_index = i * batch_size
            end_index = (i + 1) * batch_size

            # print(start_index)
            # print(end_index)

            batch_dicts = variant_update_dict[start_index:end_index]

            db.session.execute(variant_update_stmt, batch_dicts)
            db.session.commit()
    else:
        print('Update dict empty, not updating any variants...')

    # db.session.execute(variant_update_stmt, variant_update_dict) if len(variant_update_dict)>0 else print('Update dict empty, not updating any variants...')
    # db.session.commit()
    # db.session.close()