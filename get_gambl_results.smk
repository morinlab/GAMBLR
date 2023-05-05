# Sync some of the main GAMBL results used by GAMBLR to a local or cloud host for interactive R sessions etc
# This script ONLY works when you have an active VPN connection to the GSC
# It also requires an RSA key that is passphrase-protected. Other key types don't seem to work.
# You can create one and add it to your authorized_keys file if your current key is not passphase-protected

envvars:
    "GSC_USERNAME",
    "GSC_PASSPHRASE",
    "GSC_KEY"


try:
    hostname = os.environ["GSC_HOSTNAME"]
except KeyError:
    hostname = 'gphost08.bcgsc.ca'
    # use this to sync from a different host (mostly for unsynced keys or host outage)

myuser = os.environ["GSC_USERNAME"] #set the environment variable GSC_USERNAME with your GSC username
mypassphrase = os.environ["GSC_PASSPHRASE"] #set the environment variable GSC_PASSPHRASE with the passphrase for your rsa key
mySSHK = os.environ["GSC_KEY"] #  set  the environment variable GSC_KEY to the path of your RSA private key
from snakemake.remote.SFTP import RemoteProvider
SFTP = RemoteProvider(username=myuser,private_key_pass=mypassphrase,private_key=mySSHK)

configfile: "config.yml" #UPDATE: specify the full path to the GAMBLR config.yml or make a symlink so this points to it


project_base = config["default"]["project_base"]
repo_base = config["default"]["repo_base"]

flatfiles = config["default"]["results_flatfiles"]
merged = config['default']['results_merged']
versioned = config['default']['results_versioned']
resources = config["default"]["resources"]
wildcards = config["default"]["results_merged_wildcards"]

expression = merged["tidy_expression_path"]
collated = merged["collated"]

#here we specify which files are included from the GAMBLR config
db_maf = flatfiles["ssm"]["template"]["cds"]["deblacklisted"]
aug_maf = flatfiles["ssm"]["template"]["cds"]["augmented"]
full_db_maf = flatfiles["ssm"]["template"]["merged"]["deblacklisted"] + ".bgz"
full_aug_maf = flatfiles["ssm"]["template"]["merged"]["augmented"] + ".bgz"

cnv_combined = flatfiles["cnv_combined"]["template"]
sv_combined = flatfiles["sv_combined"]["template"]
blacklist = resources["blacklist"]["template"]


seq_types = list(config["default"]["unmatched_normal_ids"]["gambl"].keys())
projections = config["default"]["projections"].split(",")


rule all:
    input:
        deblacklisted = expand(db_maf,seq_type=seq_types,projection=projections),
        augmented = expand(aug_maf,seq_type=seq_types,projection=projections),
        combined_cnv = expand(cnv_combined,seq_type=['genome', 'capture'],projection=projections),
        combined_sv = expand(sv_combined,seq_type=['genome'],projection=projections),
        blacklists = expand(blacklist,seq_type=seq_types,projection=projections),
        expression = expression,
        collated = expand(collated,seq_type_filter=['genome','capture']),
        indexed_maf = expand(full_db_maf,seq_type=['genome','capture'],projection=projections),
        indexed_aug_maf = expand(full_aug_maf,seq_type=['genome','capture'],projection=projections)

#Use the relative directory for local file names (outputs) and full path for remote file names (inputs)

rule get_collated:
    input:
        collated_tsv = SFTP.remote(hostname + project_base + collated)
    output:
        collated_tsv = collated
    run:
        shell("cp {input.collated_tsv} {output.collated_tsv}")

rule get_aug_maf:
    input:
        maf = SFTP.remote(hostname + project_base + full_aug_maf),
        maf_index = SFTP.remote(hostname + project_base + full_aug_maf + ".tbi")
    output:
        maf = full_aug_maf,
        maf_index = full_aug_maf + ".tbi"
    run:
        shell("cp {input.maf_index} {output.maf_index}")
        shell("cp {input.maf} {output.maf}")

rule get_indexed_maf:
    input:
        maf = SFTP.remote(hostname + project_base + full_db_maf),
        maf_index = SFTP.remote(hostname + project_base + full_db_maf + ".tbi")
    output:
        maf = full_db_maf,
        maf_index = full_db_maf + ".tbi"
    run:
        shell("cp {input.maf_index} {output.maf_index}")
        shell("cp {input.maf} {output.maf}")

#rule get_lymphgen:
#    input:
#        lymphgen = SFTP.remote(hostname + repo_base + lymphgen)
#    output:
#        lymphgen = lymphgen
#    run:
#         shell("cp {input.lymphgen} {output.lymphgen}")

rule get_expression:
    input:
        exp = SFTP.remote(hostname + project_base + expression)
    output:
        expression = expression
    run:
        shell("cp {input.exp} {output.expression}")

rule get_mafs:
    input:
        db = SFTP.remote(hostname + project_base + db_maf),
        aug = SFTP.remote(hostname + project_base + aug_maf),
    output:
        deblacklisted = db_maf,
        db_complete = touch(db_maf + ".complete"),
        augmented = aug_maf,
        aug_complete = touch(aug_maf + ".complete")
    run:
        shell("cp {input.db} {output.deblacklisted}")
        shell("cp {input.aug} {output.augmented}")

rule get_resources:
    input:
        SFTP.remote(hostname + project_base + blacklist)
    output:
        blacklist,
        touch(blacklist + ".complete")
    run:
        shell("cp {input[0]} {output[0]}")

rule get_cnv:
    input:
        cnv = SFTP.remote(hostname + project_base + cnv_combined)
    output:
        combined_cnv = cnv_combined,
        complete = touch(cnv_combined + ".complete")
    run:
        shell("cp {input[0]} {output[0]}")

rule get_sv:
    input:
        SFTP.remote(hostname + project_base + sv_combined)
    output:
        sv_combined,
        complete = touch(sv_combined + ".complete")
    run:
        shell("cp {input[0]} {output[0]}")
