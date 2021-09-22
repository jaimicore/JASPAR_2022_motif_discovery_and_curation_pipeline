/////////////////////////////////////////////////
// Insert after the oppening of the 'body' tag // 
/////////////////////////////////////////////////

<a href="https://chart-studio.plotly.com/~jaimicore/88" target="_blank"><img class="barplot" src="barplots/TF_classes_barplot_JASPAR_2022_insects_CORE+UNVALIDATED.jpeg" ></a>
<a href="https://chart-studio.plotly.com/~jaimicore/92" target="_blank"><img class="barplot" src="barplots/TF_classes_barplot_JASPAR_2022_nematodes_CORE+UNVALIDATED.jpeg" ></a>
<a href="https://chart-studio.plotly.com/~jaimicore/90" target="_blank"><img class="barplot" src="barplots/TF_classes_barplot_JASPAR_2022_urochordata_CORE+UNVALIDATED.jpeg" ></a>
<a href="https://chart-studio.plotly.com/~jaimicore/86" target="_blank"><img class="barplot" src="barplots/TF_classes_barplot_JASPAR_2022_fungi_CORE+UNVALIDATED.jpeg" ></a>
<a href="https://chart-studio.plotly.com/~jaimicore/84" target="_blank"><img class="barplot" src="barplots/TF_classes_barplot_JASPAR_2022_plants_CORE+UNVALIDATED.jpeg" ></a>
<a href="https://chart-studio.plotly.com/~jaimicore/82" target="_blank"><img class="barplot" src="barplots/TF_classes_barplot_JASPAR_2022_vertebrates_CORE+UNVALIDATED.jpeg" ></a>





<div style='position:absolute;top:525px;float:left;padding:0px;align:center;font-size:15px;border:0px' >
<ul width='450px' style='padding-left:20px;border:0px;'>
  <li>Click on the TF name to open its corresponding link in JASPAR</li>
  <li>Click on the <i>Direct/Reverse</i> buttons to change the motif orientation</li>
  <li>Ring colors: TF classes</li>
  <li>Branch colors: Clusters assigned by <a href="https://academic.oup.com/nar/article/45/13/e119/3862068" target="_blank">RSAT matrix-clustering<a> <br>See summary</li>
  <li>Click on the barplot image above to open the interactive TF class barplot</li>
</ul>
</div>



////////////////////////////////
// Insert after image section //
////////////////////////////////
 	var innerRad_start = 200;
 	var innerRad_end = 450;
 
 	
    // Add background elements for selection
     vis.selectAll('rect_selection')
        .data(data_sample)
        .enter()
        .append('path')
        .attr('class', 'rect1')
        .attr('id', function(d){ return( 'rect_sel' + d.id_motif) })
        .attr('d', d3.svg.arc()
          .startAngle(function(d){ return( (d.start * Math.PI)/180)  }  )  //converting from degs to radians
          .endAngle(function(d){   return( (d.end   * Math.PI)/180)  }  ) //just radians
          .innerRadius(innerRad_start)         // This is the size of the donut hole
          .outerRadius(innerRad_end)
        )
        .style('fill',function(d){return(d.class)})
        .attr('transform', 'translate(0,0)')
        .style('stroke-width', '1px')
        .style('opacity', 0.15);

            
          // Color algorithm annotation layer
        vis.selectAll('annotations')
            .data(data_sample)
            .enter()
            .append('path')
            .attr('class','annotation1')
            .attr('id', function(d,i){return('path' + i) })
            .attr('d', d3.svg.arc()
              .startAngle(function(d){ return( d.start * (Math.PI/180) ) }  )  //converting from degs to radians
              .endAngle(function(d){ return( d.end * (Math.PI/180) ) }  ) //just radians
              .innerRadius(innerRad_end)         // This is the size of the donut hole
              .outerRadius(innerRad_end + (30*1))
            )
            .attr('fill', 'white')
            .attr('stroke', 'white')
            .style('stroke-width', '0px')
            .style('opacity', 1);
            

        // Color algorithm Non-Validated text layer
        vis.selectAll('annotations')
            .data(data_sample)
            .enter()
            .append('text')
            .attr('dy', 38)
            .attr('x',43)
            .style('font-size', '50px')
            .append('textPath')
            .attr('startOffset','50%')
            .style('text-anchor','middle')
            .attr('stroke','black')
            .attr('fill','black')
            .attr('xlink:href', function(d,i){return('#path' + i) })
            .text(function(d){
              var token = d.matrix_name.split('_').slice(-2).slice(0);
              token = token[0];
              if(/^UN/.test(token)){
                var text_content = d.class_nb + '*';
              } else {
                var text_content = d.class_nb;
              }
              return(text_content);
            });
            //.style('font-size', '60px');
                        
            
        // Color algorithm Non-Validated text layer
            vis.selectAll('annotations')
                .data(data_sample)
                .enter()
                .append('path')
                .attr('class','annotation2')
                .attr('d', d3.svg.arc()
                  .startAngle(function(d){ return( d.start * (Math.PI/180) ) }  )  //converting from degs to radians
                  .endAngle(function(d){ return( d.end * (Math.PI/180) ) }  ) //just radians
                  .innerRadius(innerRad_end + (30*1))         // This is the size of the donut hole
                  .outerRadius(innerRad_end + (30*2))
                )
                .attr('fill', function(d){return(d.class)})
                .attr('stroke', 'white')
                .attr('transform', 'translate(0,0)')
                .style('stroke-width', '2px')
                .style('opacity', 1);

//////////////////////////////////
// Insert before g.text section //
//////////////////////////////////



///////////////////////////////////////////////////
// Replace g.text section by the following chunk //
///////////////////////////////////////////////////
            vis.selectAll('g.text')
                .data(nodes)
                .enter()
                .append('g')
                .attr('class', 'node_text')
                .attr('transform', function(d) { return 'rotate(' + (d.x - 90) + ')translate(' + d.y + ')'; })
                .append('svg:a')
                .attr('xlink:href', function(d) { return d.link_ext; })
                .attr('target', '_blank')
                .append('text')
		.text(function(d){
              		var token = d.link_ext;
              		if(/UN/.test(token)){
                		var text_content = d.name + '*';
              		} else {
                		var text_content = d.name;
              		}
              		return d.children ? '' :  text_content;
            	})
                .attr('fill', function(d) { return d.color; })
                .style('font-size', '20px')
		.attr("font-family", "Arial")
                .attr('dx', function(d) { return d.x < 180 ? 10 : -10; })
                .attr('dy', '.31em')
                .attr('text-anchor', function(d) { return d.x < 180 ? 'start' : 'end'; })
                .attr('transform', function(d) { return d.x < 180 ? null : 'rotate(180)'; });
