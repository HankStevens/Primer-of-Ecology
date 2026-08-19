# First, unset your current PAT from .Renviron
usethis::edit_r_environ()
# In the file that opens, remove or comment out the GITHUB_PAT line
# Save and restart R
# Then set up git credentials
gitcreds::gitcreds_set()
# This will prompt you to enter a new PAT
usethis::create_github_token()
# Install development version from GitHub
devtools::install_github("mjg211/phaseR")
