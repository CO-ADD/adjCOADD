
$(document).ready(function () {
   console.log("loading htmx header")
   $('body').addEventListener("htmx:configRequest", (event) => {
       event.detail.headers["X-CSRFToken"] = "{{ csrf_token }}";
   });
});
