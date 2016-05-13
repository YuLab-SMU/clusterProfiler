# Contributing to clusterProfiler development

1. Filing a bug report or feature request in an issue.
1. Suggesting a change via a pull request.

## Issues

When filing an issue, the most important thing is to include a minimal reproducible example so that we can quickly verify the problem, and then figure out how to fix it. There are three things you need to include to make your example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, so it's easy to
    see which ones the example needs.
  
1.  The easiest way to include **data** is to use `dput()` to generate the R code 
    to recreate it. For example, to recreate the `gene` vector/list in R,
    I'd perform the following steps:
  
       1. Run `dput(gene)` in R
       2. Copy the output
       3. In my reproducible script, type `gene <- ` then paste.
       
  
1.  Spend a little bit of time ensuring that your **code** is easy for others to
    read:
  
    * make sure you've used spaces and your variable names are concise, but
      informative
  
    * use comments to indicate where your problem lies
  
    * do your best to remove everything that is not related to the problem.  
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a fresh R session and pasting your script in.


## Pull requests

