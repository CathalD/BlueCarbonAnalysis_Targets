# Add after the existing HTML content (around line 300):

# ── PATH RESOLVER ───────────────────────────────────────────────────────────────────
# Ensures working directory is BlueCarbon_Workflow_V1.0/ (project root)
# so all relative data paths (data_raw/, outputs/, etc.) resolve correctly.
local({
  # Method 1: called via source() — detect this script's location
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (!is.null(this_file)) {
    setwd(dirname(dirname(dirname(this_file))))
    return()
  }
  # Method 2: running interactively in RStudio — rstudioapi gives active file
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(active) && nchar(active) > 0) {
      setwd(dirname(dirname(dirname(normalizePath(active)))))
      return()
    }
  }
  # Method 3: fallback — warn and let the user set wd manually
  message("\u26a0 Could not auto-detect project root. Please run:\n",
          "  setwd('/path/to/BlueCarbon_Workflow_V1.0')\n",
          "before sourcing this script.")
})
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
sampling_recommendations <- ''

if (exists("sample_assessment")) {
  
  needs_more <- sample_assessment %>% filter(status != "SUFFICIENT")
  
  if (nrow(needs_more) > 0) {
    table_rows <- paste(apply(needs_more, 1, function(row) {
      priority <- ifelse(row["status"] == "INSUFFICIENT", "\U0001f534 HIGH", 
                        ifelse(row["status"] == "BELOW RECOMMENDED", "\U0001f7e1 MEDIUM", "\U0001f7e2 LOW"))
      sprintf('<tr><td>%s</td><td>%d</td><td>%s</td><td>%d</td><td>%d</td><td>%s</td></tr>',
              row["stratum"], row["n_cores"], row["status"],
              row["additional_for_minimum"], row["additional_for_recommended"], priority)
    }), collapse = "\n")
    sampling_recommendations <- paste0(
'
<h2>\u26a0\ufe0f Sampling Recommendations</h2>
<div class="highlight" style="background-color: #FFF3E0; border-left: 4px solid #F57C00;">
<p>', vm_label("<strong>Additional field sampling recommended before carbon crediting:</strong>", "<strong>Additional field sampling recommended:</strong>"), '</p>
<table>
<thead>
<tr>
<th>Stratum</th>
<th>Current Cores</th>
<th>Status</th>
<th>Additional Needed (Min)</th>
<th>Additional Needed (Recommended)</th>
<th>Priority</th>
</tr>
</thead>
<tbody>
', table_rows, '
</tbody>
</table>
<p><strong>Rationale:</strong></p>
<ul>
', vm_label("<li>VM0033 requires minimum n=3 per stratum (n=5 recommended)</li>", "<li>Statistical best practice requires minimum n=3 per stratum (n=5 recommended)</li>"), '
<li>Robust spatial models need n=10+ per stratum</li>
<li>Current sample sizes result in high uncertainty and conservative estimates</li>
<li>Additional samples will reduce uncertainty and increase creditable carbon</li>
</ul>
<p><strong>Estimated Impact:</strong></p>
<ul>
<li>Current uncertainty may reduce creditable carbon by 15-30%</li>
<li>Meeting recommended sample sizes could increase credit volume by 20-40%</li>
<li>Improved spatial coverage will strengthen verification</li>
</ul>
</div>
'
    )
  } else {
    sampling_recommendations <- paste0(
'
<h2>\u2705 Sampling Status</h2>
<div class="highlight" style="background-color: #E8F5E9; border-left: 4px solid #4CAF50;">
<p><strong>', vm_label("Sampling is adequate for carbon crediting", "Sampling is statistically adequate"), '</strong></p>
<p>', vm_label("All strata meet or exceed VM0033 recommended sample sizes (n\u22655 per stratum).", "All strata meet or exceed recommended sample sizes (n\u22655 per stratum)."), '</p>
<p>', vm_label("Proceed with verification package submission.", "Proceed with soil carbon assessment."), '</p>
</div>
'
    )
  }
}

# Insert into HTML before verification checklist section
html_content <- sprintf('
<!DOCTYPE html>
<html>
<head>
...
%s

<h2>Verification Checklist</h2>
...
',
sampling_recommendations
)