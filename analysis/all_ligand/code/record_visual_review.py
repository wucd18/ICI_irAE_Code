from common import *
pages=jr(ROOT/'review/pdf_render_checks.json');assert len(pages)==15 and all(r['render_exit_code']==0 for r in pages)
counts=pd.DataFrame(pages).groupby('file').size().to_dict();assert len(counts)==6
notes={
'Figure4_TNFSF12.pdf':'Panels a–d present, all 12 donor–family observations in each applicable donor panel, explicit atlas FDR, rendered superscript minus fixed, readable labels and legend.',
'Figure_AllLigand_Atlas.pdf':'Both alphabetical pages checked; all 79 ligand rows, all four families and three readouts shown; NE labels, separate colour scales, donor counts and FDR outlines retained; footer spacing corrected.',
'Figure_Response_Patterns.pdf':'Both pages checked; categories and counts visible; bar labels no longer touch the matrix; legend and definition text no longer overlap.',
'Figure_Reference_Benchmark.pdf':'All three pages checked; n and rho readable, every scatter retains raw axes; lowest-row labels separated from footer; TNFSF12 rings and membership labels retained.',
'Supplementary_Selected_Examples.pdf':'All five pages checked; 38 unique selected combinations balanced over pages; three actual donor markers per panel, UP/Dep/IFN direction counts in headings, no empty last-page layout.',
'Supplementary_Residual_and_AssayReference.pdf':'Both pages checked; residual correction scope separated from assay reference; donor markers, mean ticks, q values and nonsignificant open symbols visible.'}
records=[dict(file=n,pages=int(k),sha256=sha(ROOT/'figures'/n),status='PASS',review=notes[n]) for n,k in counts.items()]
jw(ROOT/'review/visual_review.json',dict(status='PASS',method='Poppler rendering of actual final PDF pages at 150 dpi, followed by model visual inspection of full-size pages and contact sheets; not merely drawing-property checks',files=records,rendered_pages=15,reviewed_utc=now()))
(ROOT/'review/VISUAL_REVIEW.md').write_text('# Final PDF visual review\n\nAll six final PDFs (15 pages) were rendered using Poppler and visually inspected.\n\n'+'\n\n'.join('## '+r['file']+'\n\n'+r['review'] for r in records)+'\n',encoding='utf8')
print('VISUAL_REVIEW_RECORDED',flush=True)
