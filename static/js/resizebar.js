$("#hidefilter").on("click", function () {
  console.log("test");
  $(this).find("i").toggleClass("bi-funnel bi-funnel-fill");
  $(".sidebar").fadeToggle("fast");
});
