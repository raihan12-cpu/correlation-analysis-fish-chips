import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
const root=path.dirname(fileURLToPath(import.meta.url));
const result=JSON.parse(await fs.readFile(path.join(root,'results/analysis_results.json'),'utf8'));
const wb=Workbook.create();
const labels={Protein_pct:'Protein (%)',Ash_pct:'Ash (%)',Fat_pct:'Fat (%)',Moisture_pct:'Moisture (%)',Linear_expansion_pct:'Expansion (%)',Reported_oil_absorption_pct:'Oil index (%)',Specific_volume_cm3_g:'Specific volume (cm³/g)',Color_L:'L*',Color_a:'a*',Color_b:'b*'};
function add(name,title,columns,rows,note){
 const s=wb.worksheets.add(name),n=columns.length;
 const col=i=>String.fromCharCode(65+i);
 s.showGridLines=false;
 s.getRangeByIndexes(0,0,rows.length+6,n).format.font={name:'Arial',size:11,color:'#202A35'};
 s.getRange('A2').values=[[title]];s.getRange('A2').format.font={size:15,bold:true};
 s.getRangeByIndexes(3,0,1,n).values=[columns];s.getRangeByIndexes(4,0,rows.length,n).values=rows;
 s.getRangeByIndexes(3,0,1,n).format={fill:'#E8EDF2',font:{bold:true},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:34,borders:{bottom:{style:'thin',color:'#637488'}}};
 s.getRangeByIndexes(4,0,rows.length,n).format.rowHeight=23;
 s.freezePanes.freezeRows(4);
 for(let c=0;c<n;c++){
  const num=rows.some(r=>typeof r[c]==='number');
  s.getRangeByIndexes(3,c,rows.length+1,1).format.columnWidth=num?16:29;
  s.getRangeByIndexes(4,c,rows.length,1).format.horizontalAlignment=num?'right':'left';
  if(num)s.getRangeByIndexes(4,c,rows.length,1).setNumberFormat(['n','Batch_number'].includes(columns[c])?'0':'0.0000');
 }
 s.getRangeByIndexes(1,n+1,1,1).values=[['Notes']];
 s.getRangeByIndexes(2,n+1,1,1).values=[[note]];
 s.getRangeByIndexes(2,n+1,1,1).format={columnWidth:75,wrapText:true,rowHeight:75};
 return s;
}
for(const [name,key,title,note] of [
 ['Pearson A-G','pearson','Pearson A–G','21 independent batches. Two-sided Pearson p-values; BH adjustment over 45 unique pairs. Pooled associations, not formulation-adjusted effects.'],
 ['Spearman A-G','spearman','Spearman A–G','21 batches. Average ranks for ties. Spearman p-values use the t approximation; interpret cautiously with small n. BH over 45 pairs separately.'],
 ['Pearson A-F','pearson_without_G','Pearson A–F','18 batches; control G excluded. Two-sided Pearson p-values; BH adjustment over 45 pairs separately.']
]){
 const rows=result[key].map(r=>[labels[r.x],labels[r.y],r.n,r.r,r.p,r.q]);
 const s=add(name,title,['Variable 1','Variable 2','n','Coefficient','p, two-sided','p, BH-adjusted'],rows,note);
 s.getRange('E5:F49').setNumberFormat('[<0.001]"<0.001";0.0000');
 for(let i=0;i<rows.length;i++)if(rows[i][5]<.05)s.getRange('F'+(i+5)).format.font={bold:true};
}
const cols=['Sample_ID','Batch_ID','Batch_number',...result.audit.variables];
add('Matched data','Matched data',cols,result.raw.map(r=>cols.map(c=>r[c])),'Source: Fish_Crackers_Correlation_PCA_HCA_Preprocessed.xlsx, Batch_data. Metadata and carbohydrate are excluded from correlation variables.');
const audit=[
 ['Independent batches',21],['Variables analyzed',10],['Unique pairs per analysis',45],['Missing values',0],
 ['Duplicate batch keys',0],['Duplicate numerical rows',result.audit.duplicate_numeric_rows],
 ['Proximate values reconciled',84],['Physical/color values reconciled',126],
 ['Carbohydrate included','No'],['PCA/HCA performed','No; only supplied preprocessing was checked'],
 ['Main Pearson BH <0.05',result.significant_counts.pearson],
 ['Pearson without G BH <0.05',result.significant_counts.without_G],
 ['Interpretation','Pooled correlations can reflect formulation differences; compositional constraints remain.'],
 ['Oil index','Weight-change-based frying index; direct oil uptake unverified.'],
 ['Heatmap scale','Fixed at -1 to +1; no significance stars.'],
 ['Probability method','Pearson: beta-distribution null. Spearman: t approximation.'],
 ['Workbook behavior','Static output; rerun Python and this builder to update.']
];
add('Audit','Data and analysis audit',['Item','Result'],audit,'No measurements were imputed or deleted. All source values retained; carbohydrate excluded from correlation calculations.');
const auditSheet=wb.worksheets.getItem('Audit');
auditSheet.getRange('B4:B21').format.columnWidth=100;
auditSheet.getRange('B5:B21').setNumberFormat('0');
auditSheet.getRange('B5:B21').format.horizontalAlignment='left';
wb.recalculate();
await fs.mkdir(path.join(root,'previews'),{recursive:true});
for(const name of ['Pearson A-G','Spearman A-G','Pearson A-F','Matched data','Audit']){
 const pic=await wb.render({sheetName:name,range:name==='Matched data'?'A1:G12':name==='Audit'?'A1:B21':'A1:F15',scale:1,format:'png'});
 await fs.writeFile(path.join(root,'previews',name.replaceAll(' ','_')+'.png'),new Uint8Array(await pic.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(root,'results/correlation_results.xlsx'));
console.log((await wb.inspect({kind:'sheet',include:'id,name'})).ndjson);

