class TreeRenderer {
  constructor(treeContainer){
    this.treeContainer = treeContainer;

    this.scale = 200;
    this.padding = 50;
    this.width = 500;
    this.height = 500;

    this.streamgraph_options = new StreamgraphOptions(100, 10);

    this.svg = this.treeContainer.append("svg")
        .attr("height", this.width)
        .attr("width", this.height)
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
    let _y = [];
    for (let i = 1; i < coords.length; i+=3) {
      _y.push(coords[i]);
    }
    const _y_max = Math.max(..._y);
    const _y_min = Math.min(..._y);

    // also change the x coordinates, so our bounding box starts at (0,0). This helps for e.g. the coordinate axis rendering
    let _x = [];
    for (let i = 0; i < coords.length; i+=3) {
      _x.push(coords[i]);
    }
    const _x_min = Math.min(..._x);

    while(coords.length > 0) {
      this.points.push(new Point(
          (coords.shift() - _x_min) * this.scale,
          ((Math.abs(_y_max - _y_min)) - (coords.shift() - _y_min)) * this.scale,
          coords.shift() * this.scale)
      );
    }

    this.x_max = Math.max(...this.points.map(p => p.x));
    this.x_min = Math.min(...this.points.map(p => p.x));
    this.y_max = Math.max(...this.points.map(p => p.y));
    this.y_min = Math.min(...this.points.map(p => p.y));

    for(let i = 0; i < this.points.length; i++) {
      this.points[i].setBranchId(branchid[i]);
      this.points[i].setKDE(kde[i]);
      this.points[i].setKDE_I0([kde_i0[i * 3], kde_i0[i * 3 + 1], kde_i0[i * 3 + 2]]);  // 3 components for each entry
      this.points[i].setKDE_I1([kde_i1[i * 3], kde_i1[i * 3 + 1], kde_i1[i * 3 + 2]]);  // 3 components for each entry
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
    this.treeRoot = this.points[this.points.findIndex(p => p.y === this.y_max)];
    this.treeRoot.drawPoint = true;

    // center the tree with viewBox
    this.svg
        .attr("viewBox", [
            this.x_min - this.padding,
            this.y_min - this.padding,
            2 * this.padding + (this.x_max - this.x_min),
            2 * this.padding + (this.y_max - this.y_min)]
        );
  }

  render(type) {
    if (!this.vtkDataSet) {
      console.error("trying to render without dataset");
      return;
    }

    this.nodelayer.empty();

    // draw the coordinate system
    this.coordinate_system
        .attr("transform", `translate(${-2 * this.padding} , 0)`);

    this.coordinate_system.append("text")
        .attr("x", -40)
        .attr("y", -20)
        .attr("font-size", 10)
        .attr("font-weight", "lighter")
        .text("Density");

    const yScale = d3.scaleLinear()
        .domain([0, 1])
        .range([this.y_max, 0]);
    const yAxis = d3.axisRight(yScale)
        .ticks(10)
        .tickSize(this.width - 2 * this.padding);
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

    // draw the streamgraph
    if (type === 'streamgraph') {
      const branches = this.points
          .map(p => p.branchId)
          .filter((v, i, self) => self.indexOf(v) === i)
          .sort((a, b) => this.xCoordOfBranchID(a) - this.xCoordOfBranchID(b))  // sort by the x coordinate of the branches
          .reverse();

      // We want each subtree to have the same order of categories (i.e. "colors") as
      // the main tree. The order is specified by the size at the root node
      // see https://stackoverflow.com/questions/3730510/javascript-sort-array-and-return-an-array-of-indices-that-indicates-the-positio
      const color_order = Array.from(Array(this.treeRoot.kde_i1.length).keys())  // create array [0, ..., this.treeRoot.kde_i1.length - 1]
          .sort((a, b) => this.treeRoot.kde_i1[b] - this.treeRoot.kde_i1[a]);  // sort the array according to the entries in this.treeRoot.kde_i1

      branches.forEach((branch, i) => {
        // calculate maximum space between the edge of current branch to the branch that is to the right of it
        const maximum_space = (i === 0)
            ? this.streamgraph_options.maxwidth_root
            : this.xCoordOfBranchID(branches[i - 1]) - this.xCoordOfBranchID(branch);

        const current_branch = this.points
            .filter(p => p.branchId === branch)
            .sort((p1, p2) => p2.y - p1.y);

        const data = current_branch
            .map(p => p.kde_i1)
            .map(kde => color_order
                .map(i => kde[i])  // reorder the kde_i1 arrays according to the indices in color_order
            )
            .map(kde => d3.cumsum(kde));

        const no_of_categories = data[0].length;

        // create a scale for mapping of the points
        const mapping = d3.scaleLinear()
            .domain([0, data[0][no_of_categories - 1]])
            .range([0, maximum_space - this.streamgraph_options.padding]);

        // draw the streamgraph
        for (let i = 0; i < no_of_categories; i++) {
          let path = d3.path();
          path.moveTo(current_branch[0].x, current_branch[0].y);
          current_branch.forEach((point, j) => {
            path.lineTo(point.x + mapping(data[j][no_of_categories - 1 - i]), point.y);
          });
          path.lineTo(current_branch[current_branch.length - 1].x, current_branch[current_branch.length - 1].y);
          path.closePath();

          this.nodelayer.append("path")
              .attr("d", path)
              .attr("fill", d3.schemeSet1[i])
              .attr("stroke", "black")
              .attr("stroke-width", 0.5);
        }
      });
    }

    // draw the lines of the graph
    const unique_branchIds = this.points
        .map(p => p.branchId)
        .filter((v, i, self) => self.indexOf(v) === i);

    unique_branchIds.forEach(bId => {
      const pointsOnBranch = this.points
          .filter(p => p.branchId === bId)
          .sort((a, b) => a.y - b.y);

      const x1 = pointsOnBranch[0].x;
      const y1 = pointsOnBranch[0].y;
      const y2 = pointsOnBranch[pointsOnBranch.length - 1].y;

      let path = d3.path();
      path.moveTo(x1, y1);
      path.lineTo(x1, y2);

      const x = this.points.indexOf(pointsOnBranch[pointsOnBranch.length - 1]);
      const _id = this.connectivityArray.filter((c, i) => i % 2 === 0).indexOf(BigInt(x));
      if (_id !== -1) {  // if _id === -1, then we are at the root node.
        const x3 = this.points[this.connectivityArray[_id * 2 + 1]].x;
        const y3 = this.points[this.connectivityArray[_id * 2 + 1]].y;
        path.lineTo(x3, y3);
      }

      this.nodelayer.append("path")
          .attr("d", path)
          .attr("fill", "none")
          .attr("stroke", "black")
          .attr("stroke-width", 2);
    });

    // draw the donut plots
    if (type === 'donut') {
      for (let i = 0; i < this.points.length; i++) {
        if (this.points[i].drawDonut) {
          this.Donut(this.points[i]);
        }
      }
    }

    // draw the points in the graph
    for (let i = 0; i < this.points.length; i++) {
      if (this.points[i].drawPoint) {
        this.nodelayer.append("circle")
            .attr("cx", this.points[i].x)
            .attr("cy", this.points[i].y)
            .attr("r", 4)
            .attr("fill", "black");
      }
    }
  }

  Donut(point) {
    // adapted from https://www.geeksforgeeks.org/d3-js-pie-function/
    let g = this.nodelayer.append("g")
        .attr("transform", `translate(${point.x}, ${point.y})`);

    const pie = d3.pie().sort(null);
    const arc = d3.arc()
        .innerRadius(8)
        .outerRadius(20);

    const arcs = g.selectAll("arcs")
        .data(pie(point.kde_i1))
        .enter()
        .append("g");

    arcs.append("path")
        .attr("fill", (data, i) => d3.schemeSet1[i])
        .attr("d", arc)
        .attr("stroke", "black")
        .style("stroke-width", 1);
  }

  xCoordOfBranchID(branchId) {
    return this.points.filter(p => p.branchId === branchId).map(p => p.x)[0];
  }
}

class Point {
  constructor(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;

    this.drawDonut = false;
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
}

class StreamgraphOptions {
  constructor(maxwidth_root, padding) {
    this.maxwidth_root = maxwidth_root;
    this.padding = padding;
  }
}
