#!/bin/bash
# vim: set ts=4 sts=4 sw=4 et:

# Crawls through the LCO archive and attempts to fit gain / readnoise from
# automatically identified pairs of sky flats and daytime biases.

# PostgreSQL Database Configuration using the LCO standard for database
# connection configuration in containerized projects.
DB_HOST="${DB_HOST:-127.0.0.1}"
DB_PORT="${DB_PORT:-5432}"
DB_NAME="${DB_NAME:-lcogt-commissioning}"
DB_USER="${DB_USER:-lcogt-commissioning}"
DB_PASS="${DB_PASS:-undefined}"
NDAYS="${NDAYS:-3}"

# SQLAlchemy database connection string
DATABASE="postgresql://${DB_USER}:${DB_PASS}@${DB_HOST}:${DB_PORT}/${DB_NAME}"

# Existing local subdirectory into which an overview HTML page
# will be rendered
OUTPUTDIR="${OUTPUTDIR:-/home/dharbeck/public_html/gainhistory}"

# Common arguments for the crawlnoisegain script
CRAWLNOISEGAIN_ARGS=(
    "--ndays=${NDAYS}"
    "--noreprocessing"
    "--loglevel=INFO"
    "--database=${DATABASE}"
)

time crawlnoisegain "${CRAWLNOISEGAIN_ARGS[@]}" --useaws --cameratype="fa" --readmode="full_frame"
time crawlnoisegain "${CRAWLNOISEGAIN_ARGS[@]}" --useaws --cameratype="fa" --readmode="central_2k_2x2"
time crawlnoisegain "${CRAWLNOISEGAIN_ARGS[@]}" --useaws --cameratype="fs" --readmode="default"
time crawlnoisegain "${CRAWLNOISEGAIN_ARGS[@]}" --useaws --cameratype="kb" --readmode="default"
time crawlnoisegain "${CRAWLNOISEGAIN_ARGS[@]}" --useaws --cameratype="ep" --readmode="MUSCAT_SLOW"
time crawlnoisegain "${CRAWLNOISEGAIN_ARGS[@]}" --useaws --cameratype="ep" --readmode="MUSCAT_FAST"


time analysegainhistory --database="${DATABASE}"

# This script doesn't check any exit codes to confirm successful operation,
# so we won't bother either, and just always assume the best and exit with
# a code of 0 (success).
exit 0
