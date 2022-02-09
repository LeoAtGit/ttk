class TreeRenderer {
  constructor(treeContainer){
    this.treeContainer = treeContainer;

    this.scale = 200;
    this.padding = 50;
    this.width = 500;
    this.height = 500;

    let topN = 3;
    this.streamgraph_options = new StreamgraphOptions(100, 10, 1, topN, "#9f9f9f");
    this.donut_options = new DonutOptions(8, 20, 10, topN, "#9f9f9f");

    this.svg = this.treeContainer.append("svg")
        .attr("height", this.height)
        .attr("width", this.width)
        .attr("style", "border:5px solid #eaeaea");

    this.root = this.svg.append("g");
    this.coordinate_system = this.root.append("g");
    this.nodelayer = this.root.append("g");

    this.vtkDataSet = null;
    this.points = [];
  }

  setVtkDataSet(dataset){
    this.vtkDataSet = dataset;

    // point information
    const _coords = this.vtkDataSet.points.coordinates.data;
    let coords = Array.from(_coords);  // convert from Float64Array to normal Array, so we can use `coords.shift()`
    const branchid = this.vtkDataSet.pointData.BranchId.data;
    const kde = this.vtkDataSet.pointData.KDE.data;
    const kde_i0 = this.vtkDataSet.pointData.KDE_ByType_I0.data;
    const kde_i1 = this.vtkDataSet.pointData.KDE_ByType_I1.data;

    // change the y coordinates, because ttk and svg use a different coordinate system
    const _y = coords.filter((d, i) => (i % 3) === 1);
    const _y_max = Math.max(..._y);
    const _y_min = Math.min(..._y);

    // also change the x coordinates, so our bounding box starts at (0,0). This helps for e.g. the coordinate axis rendering
    const _x = coords.filter((d, i) => (i % 3) === 0);
    const _x_min = Math.min(..._x);

    while(coords.length > 0) {
      this.points.push(new Point(
          (coords.shift() - _x_min) * this.scale,
          ((Math.abs(_y_max - _y_min)) - (coords.shift() - _y_min)) * this.scale,
          coords.shift() * this.scale)
      );
    }

    const nCompI0 = this.vtkDataSet.pointData.KDE_ByType_I0.nComponents;
    const nCompI1 = this.vtkDataSet.pointData.KDE_ByType_I1.nComponents;

    for(let i = 0; i < this.points.length; i++) {
      this.points[i].setBranchId(branchid[i]);
      this.points[i].setKDE(kde[i]);

      this.points[i].setKDE_I0(kde_i0.slice(i * nCompI0, (i + 1) * nCompI0));
      this.points[i].setKDE_I1(kde_i1.slice(i * nCompI1, (i + 1) * nCompI1));
    }

    // cell information
    this.connectivityArray = this.vtkDataSet.cells.connectivityArray.data;
    const offsetsArray = this.vtkDataSet.cells.offsetsArray.data;

    // sanity check that we only work with lines, and that the data for
    // the lines is 'nicely' arranged in the connectivityArray.
    //
    // If there is connectivity information for other geometric primitives,
    // this has to be implemented separately.
    for (let i = 0; i < offsetsArray.length; i++) {
      if (BigInt(i*2) !== offsetsArray[i]) {
        console.error("there were non-line segments in the data");
        console.error("alternatively, the data was not 'nicely' arranged (i.e. always i*2)");
        return;
      }
    }

    // for donut chart, we have to find out at which points we want to draw the donuts.
    // the locations are:
    //    - at the maximum of each branch
    //    - if two branches connect with each other, at the lower point we'd want a donut
    const unique_branchIds = [...new Set(this.points.map(p => p.branchId))];
    for (let i = 0; i < unique_branchIds.length; i++) {
      let branch = this.points.filter(p => p.branchId === unique_branchIds[i]);
      // point with lowest y-value is the maximum, i.e. the "highest" point in the view of that branch
      const max = Math.min(...branch.map(p => p.y));
      const pos = this.points.findIndex(p => p.branchId === unique_branchIds[i] && p.y === max);
      this.points[pos].drawDonut = true;
      this.points[pos].drawPoint = true;
    }

    for (let i = 0; i < this.connectivityArray.length; i+=2) {
      const p1 = this.points[this.connectivityArray[i]];
      const p2 = this.points[this.connectivityArray[i+1]];
      if (p1.branchId !== p2.branchId) {
        let idx = (p2.y > p1.y) ? 1 : 0;
        this.points[this.connectivityArray[i + idx]].drawDonut = true;
        this.points[this.connectivityArray[i + idx]].drawPoint = true;
      }
    }

    // at the absolute minimum (the root node of the tree), we also want to draw a point
    this.treeRoot = this.points[this.points.findIndex(p => p.y === Math.max(...this.points.map(p => p.y)))];
    this.treeRoot.drawPoint = true;
  }

  render(type) {
    this.type = type;

    if (!this.vtkDataSet) {
      console.error("trying to render without dataset");
      return;
    }

    this.nodelayer.empty();

    const tree = new Tree(this.points, this.connectivityArray, this.streamgraph_options, this.donut_options, this.nodelayer);

    // draw the streamgraph
    if (type === 'streamgraph') {
      tree.calculateLayoutStreamgraph();
      tree.moveToOrigin();
      tree.rescale(this.width, this.height, type);
      this.centerView(tree);
      this.drawCoordinateSystem();
      tree.drawStreamGraph();

      tree.drawEdges();
      tree.drawPoints();
    }

    // draw the donut plots
    if (type === 'donut') {
      tree.calculateLayoutDonut();
      tree.moveToOrigin();
      tree.rescale(this.width, this.height, type);
      this.centerView(tree);
      this.drawCoordinateSystem();

      tree.drawEdges();
      tree.drawDonut();
      tree.drawPoints();
    }
  }

  centerView(tree) {
    // change x_layout and y_layout coordinates such that they are centered for this particular viewbox configuration
    const x_max = Math.max(...this.points.map(p => p.x_layout))
        + ((this.type === 'streamgraph') ? this.streamgraph_options.maxwidth_root : this.donut_options.outerRadius);
    const y_max = Math.max(...this.points.map(p => p.y_layout));

    if (x_max >= this.width) {
      tree.all_points.forEach(p => p.y_layout += (this.height - y_max) / 2)
    }

    if (y_max >= this.height) {
      tree.branches.forEach(b => b.setXLayout(b.x_layout + (this.width - x_max) / 2));
    }

    this.svg
        .attr("viewBox", [
          -this.padding,
          -this.padding,
          this.width + 2 * this.padding,
          this.height + 2 * this.padding
          ]
        );
  }

  drawCoordinateSystem() {
    this.coordinate_system.append("text")
        .attr("x", -30)
        .attr("y", this.height / 2)
        .attr("font-size", 10)
        .attr("font-weight", "lighter")
        .attr("transform", `rotate(-90, ${-30}, ${this.height / 2})`)
        .text("Density");

    const y_max = Math.max(...this.points.map(p => p.y_layout));
    const y_min = Math.min(...this.points.map(p => p.y_layout));

    const yScale = d3.scaleLinear()
        .domain([0, 1])
        .range([y_max, y_min]);
    const yAxis = d3.axisRight(yScale)
        .ticks(10)
        .tickSize(this.width + this.padding);
    const gY = this.coordinate_system.append("g")
        .call(yAxis)
        .call(g => g.select(".domain")
            .remove())
        .call(g => g.selectAll(".tick line")
            .attr("stroke-opacity", 0.25))
        .call(g => g.selectAll(".tick text")
            .attr("x", -20))
        .attr("font-weight", "lighter");

    this.svg.call(d3.zoom()
        // .extent([[0, 0], [500, 500]])
        // .scaleExtent([1, 8])
        .on("zoom", ({transform}) => {
          this.nodelayer.attr("transform", transform);
          gY.call(yAxis.scale(transform.rescaleY(yScale)))
              .call(g => g.select(".domain")
                  .remove())
              .call(g => g.selectAll(".tick line")
                  .attr("stroke-opacity", 0.25))
              .call(g => g.selectAll(".tick text")
                  .attr("x", -20))
              .attr("font-weight", "lighter");
        })
    );
  }
}

function getSortedIndices(arr) {
  // see https://stackoverflow.com/questions/3730510/javascript-sort-array-and-return-an-array-of-indices-that-indicates-the-positio
  return Array.from(Array(arr.length).keys())  // create array [0, ..., arr.length - 1]
      .sort((a, b) => arr[b] - arr[a]);  // sort the array according to the entries in arr
}

function sum(arr) {
  return arr.reduce((c, k) => c + k, 0);
}

class Tree {
  constructor(points, connectivity, streamgraph_options, donut_options, nodelayer) {
    this.all_points = points;
    this.connectivity = connectivity;
    this.streamgraph_options = streamgraph_options;
    this.donut_options = donut_options;
    this.nodelayer = nodelayer;

    this.no_of_categories = this.all_points[0].kde_i1.length;

    this.unique_branchIds = this.all_points
        .map(p => p.branchId)
        .filter((v, i, self) => self.indexOf(v) === i)
        .sort((a, b) => this.xCoordFromBranchId(a) - this.xCoordFromBranchId(b))  // sort by the x coordinate of the branches
        .reverse();

    this.treeRoot = this.all_points[this.all_points.findIndex(p => p.y === Math.max(...this.all_points.map(p => p.y)))];

    // We want each subtree to have the same order of categories (i.e. "colors") as
    // the main tree. The order is specified by the size at the root node
    this.color_order = getSortedIndices(this.treeRoot.kde_i1);

    this.all_points.forEach(p => p.reorderKDE_I1(this.color_order, this.streamgraph_options.topN));

    // create a scale for mapping of the points, which is shared at each subbranch
    this.mapping = d3.scaleLinear()
        .domain([0, sum(this.treeRoot.kde_i1)])
        .range([0, this.streamgraph_options.maxwidth_root - this.streamgraph_options.padding]);

    this.branches = this.unique_branchIds
        .map((id, i) => new Branch(this, id, i === 0));
  }

  calculateLayoutStreamgraph() {
    this.branches.forEach(b => {
      if (!b.is_root_branch) {
        const space_needed = this.mapping(b.bottom.kde_i1_sorted_and_ordered_cumsum[this.no_of_categories - 1]) + this.streamgraph_options.padding;
        const new_x = this.findNeighborBranch(b).x_layout - space_needed;
        b.setXLayout(new_x);
      }
    });
  }

  calculateLayoutDonut() {
    this.branches.forEach(b => {
      if (!b.is_root_branch) {
        const space_needed = this.donut_options.outerRadius * 2 + this.donut_options.padding;
        const new_x = this.findNeighborBranch(b).x_layout - space_needed;
        b.setXLayout(new_x);
      }
    });
  }

  moveToOrigin() {
    const _x_min = Math.min(...this.all_points.map(p => p.x_layout));
    this.branches.forEach(b => b.setXLayout(b.x_layout - _x_min));
  }

  rescale(width, height, type) {
    const _x_max = Math.max(...this.all_points.map(p => p.x_layout))
        + ((type === 'streamgraph') ? this.streamgraph_options.maxwidth_root : this.donut_options.outerRadius);
    const _y_max = Math.max(...this.all_points.map(p => p.y_layout));

    // calculate the two possible scaling factors
    let f1 = width / _x_max;
    let f2 = height / _y_max;

    // check which scaling factors keep us in the bounds of width * height
    if (f1 * _y_max > height)
      f1 = 0;

    if (f2 * _x_max > width)
      f2 = 0;

    if (f1 === 0 && f2 === 0) {
      console.error("something went wrong trying to scale the values. This should never happen I think...")
      return;
    }

    let scaling_factor = 0;
    // check with which scaling factor we would get the largest bounding box
    if (width * _y_max * f1 > height * _x_max * f2) {
      scaling_factor = f1;
    } else {
      scaling_factor = f2;
    }

    if (scaling_factor === 0) {
      console.error("something went wrong trying to scale the values. This should never happen I think...")
      return;
    }

    // actually do the rescaling
    this.branches.forEach(b => b.setXLayout(b.x_layout * scaling_factor));
    this.all_points.forEach(p => p.y_layout *= scaling_factor);

    // we also have to apply the scaling to the mapping!
    this.mapping.range([0, (this.streamgraph_options.maxwidth_root - this.streamgraph_options.padding) * scaling_factor]);
    this.streamgraph_options.maxwidth_root *= scaling_factor;
  }

  drawStreamGraph() {
    this.branches.forEach(b => b.drawStreamGraph());
  }

  drawDonut() {
    // adapted from https://www.geeksforgeeks.org/d3-js-pie-function/
    this.all_points.filter(p => p.drawDonut).forEach(p => {
      let g = this.nodelayer.append("g")
          .attr("transform", `translate(${p.x_layout}, ${p.y_layout})`);

      const pie = d3.pie();
      const arc = d3.arc()
          .innerRadius(this.donut_options.innerRadius)
          .outerRadius(this.donut_options.outerRadius);

      let pie_data = p.kde_i1_sorted.slice(0, this.donut_options.topN);
      pie_data.push(sum(p.kde_i1_sorted.slice(this.donut_options.topN)));

      const arcs = g.selectAll("arcs")
          .data(pie(pie_data))
          .enter()
          .append("g");

      arcs.append("path")
          .attr("fill", (data, i) => {
                if (i !== this.donut_options.topN)
                  return d3.schemeSet1[this.color_order.indexOf(p.kde_i1_sorted_indices[i])];
                else
                  return this.donut_options.color_of_ignored;
              })
          .attr("d", arc)
          .attr("stroke", "black")
          .style("stroke-width", 1);
    });
  }

  drawEdges() {
    this.branches.forEach(b => {
      let path = d3.path();

      path.moveTo(b.x_layout, b.top.y_layout);
      path.lineTo(b.x_layout, b.connecting_point.y_layout);  // normally to b.bottom.y_layout, but we want to trick the layout drawing a bit
      path.lineTo(b.connecting_point.x_layout, b.connecting_point.y_layout);

      this.nodelayer.append("path")
          .attr("d", path)
          .attr("fill", "none")
          .attr("stroke", "black")
          .attr("stroke-width", this.streamgraph_options.edgeWidth);
    });
  }

  drawPoints() {
    this.all_points.filter(p => p.drawPoint).forEach(p =>
      this.nodelayer.append("circle")
          .attr("cx", p.x_layout)
          .attr("cy", p.y_layout)
          .attr("r", 4)
          .attr("fill", "black")
    );
  }

  findNeighborBranch(branch) {
    // TODO this might need to be more sophisticated, if there are branches with the same x-coordinates
    return this.branches[this.branches.indexOf(branch) - 1];
  }

  xCoordFromBranchId(bId) {
    return this.all_points.filter(p => p.branchId === bId).map(p => p.x)[0];
  }

  xLayoutCoordFromBranchId(bId) {
    return this.all_points.filter(p => p.branchId === bId).map(p => p.x_layout)[0];
  }
}

class Branch {
  constructor(tree, branchId, is_root_branch) {
    this.tree = tree;
    this.branchId = branchId;
    this.is_root_branch = is_root_branch;  // boolean

    this.x = this.tree.xCoordFromBranchId(this.branchId);
    this.x_layout = this.tree.xLayoutCoordFromBranchId(this.branchId);

    this.branch_points_sorted = this.tree.all_points
        .filter(p => p.branchId === this.branchId)
        .sort((a, b) => a.y - b.y);

    this.top = this.branch_points_sorted[0];
    this.bottom = this.branch_points_sorted[this.branch_points_sorted.length - 1];

    if (!this.is_root_branch) {
      const __tmp = this.tree.all_points.indexOf(this.bottom);
      const __id = this.tree.connectivity.filter((c, i) => i % 2 === 0).indexOf(BigInt(__tmp));
      // the point to which the branch is connected to
      this.connecting_point = this.tree.all_points[this.tree.connectivity[__id * 2 + 1]];
    } else {
      this.connecting_point = this.bottom;
    }
  }

  setXLayout(x_layout) {
    this.x_layout = x_layout;
    this.branch_points_sorted.forEach(p => p.setXLayout(x_layout));
  }

  drawStreamGraph() {
    for (let i = this.tree.no_of_categories - 1; i >= 0; i--) {
      let path = d3.path();
      path.moveTo(
          this.top.x_layout + this.tree.streamgraph_options.edgeWidth / 2,
          this.top.y_layout
      );
      this.branch_points_sorted.forEach(point => {
        path.lineTo(
            point.x_layout + this.tree.mapping(point.kde_i1_sorted_and_ordered_cumsum[i])
              + this.tree.streamgraph_options.edgeWidth / 2,
            point.y_layout
        );
      });
      const __lastPoint = this.branch_points_sorted[this.branch_points_sorted.length - 1];
      path.lineTo(
          __lastPoint.x_layout
            + this.tree.mapping(__lastPoint.kde_i1_sorted_and_ordered_cumsum[i])
            + this.tree.streamgraph_options.edgeWidth / 2,
          this.connecting_point.y_layout
      );
      path.lineTo(
          this.bottom.x_layout + this.tree.streamgraph_options.edgeWidth / 2,
          this.connecting_point.y_layout
      );
      path.closePath();

      let color = this.tree.streamgraph_options.color_of_ignored;
      if (i < this.tree.streamgraph_options.topN) {
        color = d3.schemeSet1[this.tree.color_order.indexOf(this.bottom.topN_indices_ordered[i]) % 9];
      }

      this.tree.nodelayer.append("path")
          .attr("d", path)
          .attr("fill", color);
    }
  }
}

class Point {
  constructor(x, y, z) {
    this.x = x;
    this.x_layout = x;
    this.y = y;
    this.y_layout = y;
    this.z = z;

    this.drawDonut = false;
  }

  reorderKDE_I1(order, topN) {
    // sort kde_i1 by size
    this.kde_i1_sorted_indices = getSortedIndices(this.kde_i1);
    this.kde_i1_sorted = this.kde_i1_sorted_indices.map(i => this.kde_i1[i]);

    // get the indices of the largest topN values
    this.topN_indices = this.kde_i1_sorted_indices.slice(0, topN);
    this.other_indices = this.kde_i1_sorted_indices.slice(topN);
    // reorder those indices according to the `order` array
    this.topN_indices_ordered = this.topN_indices.sort((i, j) => order.indexOf(i) - order.indexOf(j));

    this.kde_i1_sorted_and_ordered = this.topN_indices_ordered.map(i => this.kde_i1[i]);
    this.kde_i1_sorted_and_ordered = this.kde_i1_sorted_and_ordered.concat(this.other_indices.map(i => this.kde_i1[i]));

    this.kde_i1_sorted_and_ordered_cumsum = d3.cumsum(this.kde_i1_sorted_and_ordered);
  }

  setBranchId(id) {
    this.branchId = id;
  }

  setKDE(kde) {
    this.kde = kde;
  }

  setKDE_I0(kde) {
    this.kde_i0 = kde;
  }

  setKDE_I1(kde) {
    this.kde_i1 = kde;
  }

  setXLayout(x) {
    this.x_layout = x;
  }

  setYLayout(y) {
    this.y_layout = y;
  }
}

class StreamgraphOptions {
  constructor(maxwidth_root, padding, edgeWidth, topN, color_of_ignored) {
    this.maxwidth_root = maxwidth_root;
    this.padding = padding;
    this.edgeWidth = edgeWidth;
    this.topN = topN;
    this.color_of_ignored = color_of_ignored;
  }
}

class DonutOptions {
  constructor(innerRadius, outerRadius, padding, topN, color_of_ignored) {
    this.innerRadius = innerRadius;
    this.outerRadius = outerRadius;
    this.padding = padding;
    this.topN = topN;
    this.color_of_ignored = color_of_ignored;
  }
}
