class GraphOptions {
  constructor() {
  }

  setDonutGraph(donut) {
    this.donutGraph = donut;
  }

  setPieGraph(pie) {
    this.pieGraph = pie;
  }

  setStreamGraphAbs(stream) {
    this.streamGraphAbs = stream;
  }

  setStreamGraphRel(stream) {
    this.streamGraphRel = stream;
  }

  setStreamGraphNoScale(stream) {
    this.streamGraphNoScale = stream;
  }

  render() {
    d3.select("#color_scheme").on("input", e => {
      this.streamGraphAbs.changeColorScheme(e.target.value);
      this.streamGraphAbs.render("streamgraph");
      this.streamGraphAbs.doTransform();

      this.streamGraphRel.changeColorScheme(e.target.value);
      this.streamGraphRel.render("streamgraph");
      this.streamGraphRel.doTransform();

      this.streamGraphNoScale.changeColorScheme(e.target.value);
      this.streamGraphNoScale.render("streamgraph");
      this.streamGraphNoScale.doTransform();

      this.donutGraph.changeColorScheme(e.target.value);
      this.donutGraph.render("donut");
      this.donutGraph.doTransform();

      this.pieGraph.changeColorScheme(e.target.value);
      this.pieGraph.render("pie");
      this.pieGraph.doTransform();

      // so the color legend will be updated as well
      d3.select("#topN").node().dispatchEvent(new Event("input"));
    });

    d3.select("#topN").on("input", e => {
      const topN = e.target.value;
      d3.select("#label_topN").text(`${topN}`);

      this.streamGraphAbs.streamgraph_options.topN = parseInt(topN);
      this.streamGraphAbs.render("streamgraph");
      this.streamGraphAbs.doTransform();

      this.streamGraphRel.streamgraph_options.topN = parseInt(topN);
      this.streamGraphRel.render("streamgraph");
      this.streamGraphRel.doTransform();

      this.streamGraphNoScale.streamgraph_options.topN = parseInt(topN);
      this.streamGraphNoScale.render("streamgraph");
      this.streamGraphNoScale.doTransform();

      this.donutGraph.donut_options.topN = parseInt(topN);
      this.donutGraph.render("donut");
      this.donutGraph.doTransform();

      this.pieGraph.donut_options.topN = parseInt(topN);
      this.pieGraph.render("pie");
      this.pieGraph.doTransform();

      // render the color legend
      let labels = this.donutGraph.tree.getLabels();

      // transform the labels to an array that is sorted by color_order
      let labels_array = [];
      for (let key in labels) {
        labels_array.push([key, labels[key]])
      }
      labels_array = labels_array.sort((a, b) =>
        this.donutGraph.tree.color_order.indexOf(parseInt(a[0])) - this.donutGraph.tree.color_order.indexOf(parseInt(b[0])));

      d3.selectAll("#color_legend .entry").remove();

      labels_array.forEach(d => {
        let entry = d3.select("#color_legend")
          .append("div")
          .attr("class", "entry row")

        entry.append("div")
          .attr("style", `background-color: ${d[1]}; border: solid #aaaaaa`)
          .attr("class", "col-1 rounded-circle");

        entry.append("div")
          .attr("class", "col-10")
          .text(this.donutGraph.labels[d[0]]);
      });

      // add the "Other" label to the legend
      let entry = d3.select("#color_legend")
        .append("div")
        .attr("class", "entry row")

      entry.append("div")
        .attr("style", `background-color: ${this.donutGraph.donut_options.color_of_ignored}; border: solid #aaaaaa`)
        .attr("class", "col-1 rounded-circle");

      entry.append("div")
        .attr("class", "col-10")
        .text("Other");
    });

    d3.select("#topN_map").on("input", e => {
      let topN_map = parseInt(e.target.value);
      d3.select("#label_topN_map").text(`${topN_map}`);

      topN_map -= 1;

      this.streamGraphAbs.render("streamgraph", topN_map);
      this.streamGraphAbs.doTransform();
      this.streamGraphRel.render("streamgraph", topN_map);
      this.streamGraphRel.doTransform();
      this.streamGraphNoScale.render("streamgraph", topN_map);
      this.streamGraphNoScale.doTransform();
      this.donutGraph.render("donut", topN_map);
      this.donutGraph.doTransform();
      this.pieGraph.render("pie", topN_map);
      this.pieGraph.doTransform();
      kdeRenderer.computeMaskNoSelection();
      kdeRenderer.update_render();
    });

    d3.select("#graphScale").on("input", e => {
      const maxwidthRoot = e.target.value;
      d3.select("#label_graphScale").text(`${maxwidthRoot}%`);

      this.streamGraphAbs.streamgraph_options.maxwidth_root
        = this.streamGraphAbs.streamgraph_options.maxwidth_root_max * (parseInt(maxwidthRoot) / 100);
      this.streamGraphAbs.resetLayoutingCoords();
      this.streamGraphAbs.render("streamgraph");
      this.streamGraphAbs.doTransform();
    });

    d3.select("#treeScale").on("input", e => {
      const maxwidthRoot = e.target.value;
      d3.select("#label_treeScale").text(`${maxwidthRoot}`);

      this.streamGraphAbs.streamgraph_options.treeScale = parseInt(maxwidthRoot);
      this.streamGraphAbs.resetLayoutingCoords();
      this.streamGraphAbs.render("streamgraph");
      this.streamGraphAbs.doTransform();

      this.streamGraphRel.streamgraph_options.treeScale = parseInt(maxwidthRoot);
      this.streamGraphRel.resetLayoutingCoords();
      this.streamGraphRel.render("streamgraph");
      this.streamGraphRel.doTransform();

      this.streamGraphNoScale.streamgraph_options.treeScale = parseInt(maxwidthRoot);
      this.streamGraphNoScale.resetLayoutingCoords();
      this.streamGraphNoScale.render("streamgraph");
      this.streamGraphNoScale.doTransform();

      this.donutGraph.streamgraph_options.treeScale = parseInt(maxwidthRoot);
      this.donutGraph.resetLayoutingCoords();
      this.donutGraph.render("donut");
      this.donutGraph.doTransform();

      this.pieGraph.streamgraph_options.treeScale = parseInt(maxwidthRoot);
      this.pieGraph.resetLayoutingCoords();
      this.pieGraph.render("pie");
      this.pieGraph.doTransform();
    });

    d3.select("#contour_count").on("input", e => {
      kdeRenderer.update_nContours(parseInt(e.target.value));
      d3.select("#label_contour_count").text(`${e.target.value}`);
    });

    d3.select("#hatching_count").on("input", e => {
      kdeRenderer.update_nHatching(parseInt(e.target.value));
      d3.select("#label_hatching_count").text(`${e.target.value}`);
    });

    d3.select("#hatching_width").on("input", e => {
      kdeRenderer.update_HatchingWidth(parseInt(e.target.value));
      d3.select("#label_hatching_width").text(`${e.target.value}`);
    });

    d3.select("#graph_type").on("input", e => {
      if (e.target.value === "streamgraph_abs") {
        this.streamGraphRel.doTransform();
        this.streamGraphNoScale.doTransform();
        this.streamGraphAbs.doTransform();
        this.donutGraph.doTransform();
        this.pieGraph.doTransform();
        d3.select("#treeContainerDonut").attr("style", "display: none");
        d3.select("#treeContainerPie").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphAbs").attr("style", "");
        d3.select("#treeContainerStreamgraphRel").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphNoScale").attr("style", "display: none");
      }

      if (e.target.value === "streamgraph_rel") {
        this.streamGraphRel.doTransform();
        this.streamGraphNoScale.doTransform();
        this.streamGraphAbs.doTransform();
        this.donutGraph.doTransform();
        this.pieGraph.doTransform();
        d3.select("#treeContainerDonut").attr("style", "display: none");
        d3.select("#treeContainerPie").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphAbs").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphNoScale").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphRel").attr("style", "");
      }

      if (e.target.value === "streamgraph_no_scale") {
        this.streamGraphRel.doTransform();
        this.streamGraphNoScale.doTransform();
        this.streamGraphAbs.doTransform();
        this.donutGraph.doTransform();
        this.pieGraph.doTransform();
        d3.select("#treeContainerDonut").attr("style", "display: none");
        d3.select("#treeContainerPie").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphAbs").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphNoScale").attr("style", "");
        d3.select("#treeContainerStreamgraphRel").attr("style", "display: none");
      }

      if (e.target.value === "donut") {
        this.streamGraphRel.doTransform();
        this.streamGraphNoScale.doTransform();
        this.streamGraphAbs.doTransform();
        this.donutGraph.doTransform();
        this.pieGraph.doTransform();
        d3.select("#treeContainerDonut").attr("style", "");
        d3.select("#treeContainerPie").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphAbs").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphRel").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphNoScale").attr("style", "display: none");
      }

      if (e.target.value === "pie") {
        this.streamGraphRel.doTransform();
        this.streamGraphNoScale.doTransform();
        this.streamGraphAbs.doTransform();
        this.donutGraph.doTransform();
        this.pieGraph.doTransform();
        d3.select("#treeContainerDonut").attr("style", "display: none");
        d3.select("#treeContainerPie").attr("style", "");
        d3.select("#treeContainerStreamgraphAbs").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphRel").attr("style", "display: none");
        d3.select("#treeContainerStreamgraphNoScale").attr("style", "display: none");
      }
    });

    d3.select("#more_toggle").on("click", e => {
      if (d3.selectAll(".more").attr("style") === "display: none") {
        d3.selectAll(".more").attr("style", "display: inherit");
      } else {
        d3.selectAll(".more").attr("style", "display: none");
      }
    });

    // so the values are displayed when the side is loaded
    d3.select("#graphScale").node().dispatchEvent(new Event("input"));
    d3.select("#treeScale").node().dispatchEvent(new Event("input"));
    d3.select("#topN").node().dispatchEvent(new Event("input"));
    d3.select("#color_scheme").node().dispatchEvent(new Event("input"));
    d3.select("#contour_count").node().dispatchEvent(new Event("input"));
    d3.select("#hatching_count").node().dispatchEvent(new Event("input"));
    d3.select("#hatching_width").node().dispatchEvent(new Event("input"));
    d3.select("#graph_type").node().dispatchEvent(new Event("input"));
    d3.select("#topN_map").node().dispatchEvent(new Event("input"));
  }
}