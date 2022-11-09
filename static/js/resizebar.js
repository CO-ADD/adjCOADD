const resizer = document.querySelector("#resizer");
const sidebar = document.querySelector("#sidebar");

resizer.addEventListener("mousedown", (event) => {
    document.addEventListener("mousemove", resize, false);
    document.addEventListener(
        "mouseup",
        () => {
            document.removeEventListener("mousemove", resize, false);
        },
        false
    );
});

function resize(e) {
    var sdsize = Math.min(500, e.x);
    var size = `${sdsize}px`;
    sidebar.style.flexBasis = size;
}