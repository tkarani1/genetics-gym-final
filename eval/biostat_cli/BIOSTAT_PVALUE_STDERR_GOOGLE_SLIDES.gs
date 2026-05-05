/**
 * Creates a Google Slides presentation: biostat_cli p-values, stderr, bootstrap.
 *
 * HOW TO RUN
 * 1. Open https://script.google.com → New project.
 * 2. Paste this entire file as Code.gs (replace default content).
 * 3. Save. Select function `createBiostatPresentation` and click Run.
 * 4. Authorize (Slides + Drive scope). Check Execution log for the presentation URL.
 */

function createBiostatPresentation() {
  var pres = SlidesApp.create('biostat_cli — P-values, stderr & bootstrap');
  // Remove the default first slide so we control layout uniformly.
  pres.getSlides()[0].remove();

  var slides = [
    {
      title: 'biostat_cli',
      body:
        'P-values, standard errors & bootstrap\n\n' +
        'Focus: variant-level enrichment & rate ratio at percentile thresholds\n\n' +
        'Code: biostat_cli/stats/binary.py, biostat_cli/cli.py',
    },
    {
      title: 'Contingency at threshold t',
      body:
        'Case = is_pos true,  Control = is_pos false\n' +
        'Above = score >= t  (t in [0,1] percentile fraction)\n\n' +
        '2×2 table (rows: above / below t, columns: case / control):\n\n' +
        '        Case    Control\n' +
        'Above   TP      FP\n' +
        'Below   FN      TN',
    },
    {
      title: 'Point estimates',
      body:
        'Enrichment:\n' +
        '  Enr = [TP/(TP+FN)] / [FP/(FP+TN)]\n\n' +
        'Rate ratio (per-eval cohort totals N1, N2 from CLI or resources JSON):\n' +
        '  RR = [TP/N1] / [FP/N2]',
    },
    {
      title: 'Fisher — p-value',
      body:
        'Fisher exact test on the 2×2 table, one-sided alternative “greater”.\n\n' +
        'Same contingency drives enrichment and rate_ratio p-values in this mode.\n\n' +
        'Interpretation: evidence that cases are over-represented among scores above t vs independence.',
    },
    {
      title: 'Fisher — std_error column',
      body:
        'Enrichment:\n' +
        '  std_error is NOT the SE of Enr on the ratio scale.\n' +
        '  It is the Wald SE of log(OR) with Haldane–Anscombe +0.5 on each cell:\n' +
        '  SE_logOR = sqrt( 1/(TP+0.5) + 1/(FP+0.5) + 1/(FN+0.5) + 1/(TN+0.5) )\n\n' +
        'Rate ratio:\n' +
        '  std_error = NaN when p-value method is Fisher.',
    },
    {
      title: 'Poisson — p-value',
      body:
        'Among rows below t:  pi_below = FN / (FN + TN)\n\n' +
        'Let m = TP + FP (rows above t).\n' +
        'Under the heuristic null: expected cases above t is\n' +
        '  lambda = pi_below * m\n\n' +
        'One-sided upper tail vs Poisson(lambda):\n' +
        '  p = P(X >= TP),  X ~ Poisson(lambda)\n' +
        '(implemented as Poisson survival from TP−1).',
    },
    {
      title: 'Poisson — std_error column',
      body:
        'Enrichment: std_error = NaN (no analytic SE for Enr in this mode).\n\n' +
        'Rate ratio: delta method on log RR\n' +
        '  SE_logRR = sqrt(1/TP + 1/FP)\n' +
        '  SE_RR ≈ RR * SE_logRR\n\n' +
        'Requires TP > 0, FP > 0, and finite RR.',
    },
    {
      title: 'Bootstrap (--bootstrap N)',
      body:
        'Resample rows with replacement N times; recompute stats each time.\n\n' +
        'Unchanged on full data:\n' +
        '  value, p_value  (still Fisher or Poisson from the original table)\n\n' +
        'Overwritten:\n' +
        '  std_error ← sample SD of replicate “value” across bootstraps (ddof=1)\n\n' +
        'If N < 2 or too few finite replicates → std_error = NaN.',
    },
    {
      title: 'Bootstrap — formula',
      body:
        'For each output row (same eval, filter, score, threshold, stat):\n\n' +
        '  theta_bar = (1/N) * sum_b theta_hat_b\n' +
        '  s^2 = (1/(N-1)) * sum_b (theta_hat_b - theta_bar)^2\n\n' +
        'std_error = s',
    },
    {
      title: 'Comparison',
      body:
        '                Fisher              Poisson\n' +
        'p-value         Exact Fisher        Poisson tail vs lambda\n' +
        'Enr std_error   SE(log OR), +0.5    NaN\n' +
        'RR std_error    NaN                 RR * sqrt(1/TP + 1/FP)\n\n' +
        'Bootstrap: empirical SD of value; does NOT replace p-value.',
    },
    {
      title: 'CLI',
      body:
        'biostat_cli:\n' +
        '  --pvalue-method fisher | poisson\n' +
        '  --bootstrap [N]  (default N=100 if flag alone)\n\n' +
        'figure1-pipeline passes --pvalue-method and --bootstrap into the same run logic.',
    },
    {
      title: 'Caveats',
      body:
        'Fisher vs Poisson changes the p-value definition, not the Enr/RR formulas.\n\n' +
        'Enrichment Fisher std_error is on log-odds scale, not Enr itself.\n\n' +
        'Bootstrap reflects row resampling of the prepared frame (no clustering unless encoded in the table).',
    },
    {
      title: 'Questions',
      body: 'See: genetics-gym-final/eval/biostat_cli/stats/binary.py',
    },
  ];

  for (var i = 0; i < slides.length; i++) {
    appendTitleBodySlide(pres, slides[i].title, slides[i].body);
  }

  var url = pres.getUrl();
  Logger.log('Presentation URL: ' + url);
  // Browser alert only works in some contexts; log is reliable.
  return url;
}

/**
 * @param {GoogleAppsScript.Slides.Presentation} pres
 * @param {string} title
 * @param {string} body
 */
function appendTitleBodySlide(pres, title, body) {
  var slide = pres.appendSlide(SlidesApp.PredefinedLayout.TITLE_AND_BODY);
  slide.getPlaceholder(SlidesApp.PlaceholderType.TITLE).asShape().getText().setText(title);
  slide.getPlaceholder(SlidesApp.PlaceholderType.BODY).asShape().getText().setText(body);
}
