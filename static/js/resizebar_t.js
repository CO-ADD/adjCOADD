var resizer = [...document.querySelectorAll(".resizer")];
var sidebar = [...document.querySelectorAll(".sidebar")];

resizer.forEach((i) => i.addEventListener("mousedown", (event) => {
    console.log("enve")
    document.addEventListener("mousemove", resize, false);
    document.addEventListener(
        "mouseup",
        () => {
            document.removeEventListener("mousemove", resize, false);
        },
        false
    );
}));

function resize(e) {
    var sdsize = Math.min(500, e.x);
    var size = `${sdsize}px`;
    sidebar.forEach((i) => i.style.flexBasis = size);
}

