class GraphOptions {
  constructor() {
  }

  setDonutGraph(donut) {
    this.donutGraph = donut;
  }

  setStreamGraph(stream) {
    this.streamGraph = stream;
  }

  render() {
    // TopN
    d3.select("#topN").on("input", d => {
      const topN = d.target.value;
      d3.select("#label_topN").text(`${topN} TopN Categories`);

      this.streamGraph.streamgraph_options.topN = parseInt(topN);
      this.streamGraph.render("streamgraph");
      this.streamGraph.doTransform();

      this.donutGraph.donut_options.topN = parseInt(topN);
      this.donutGraph.render("donut");
      this.donutGraph.doTransform();

      // get all of the color mapping
      console.log(this.donutGraph.tree.getLabels());
      console.log(this.donutGraph.tree.color_mapping);
      console.log(this.donutGraph.labels);
    });
    d3.select("#label_topN").text(`${d3.select("#topN").node().value} TopN Categories`);

    // d3.select("#color_legend").append()
  }
}