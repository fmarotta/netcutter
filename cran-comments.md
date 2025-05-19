## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.

* Checks complain about a possibly invalid URL to stackoverflow, possibly because the site is using Cloudflare to verify that a human is accessing it.

* Skipping checking HTML validation: no command 'tidy' found

## Addressing CRAN comments

* References in the DESCRIPTION were formatted as requested

* Redundant R in the DESCRIPTION was removed

* T and F were converted to TRUE and FALSE

* The global environment is no longer modified
