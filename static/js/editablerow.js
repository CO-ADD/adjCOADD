$(document).ready(function () {
  $(document).on("dblclick", ".editablerow", function () {
    console.log("editable...");
    var value = $(this).text();
    var input =
      "<input type='text' class='input-data' value='" + value + "' />";
    $(this).html(input);
    $(this).removeClass("editablerow");
  });

  $(document).on("blur", ".input-data", function () {
    var value = $(this).val();
    var td = $(this).parent("td");
    $(this).remove();
    td.html(value);
    td.addClass("editablerow");
    var name = td.data("name");
    var type = td.data("type");
    data = { dict_value: name, type: type, value: value };
    console.log(data);
    sendToServer(data);
  });
  $(document).on("keypress", ".input-data", function (e) {
    if (e.keyCode === 13) {
      var value = $(this).val();
      var td = $(this).parent("td");
      $(this).remove();
      td.html(value);
      td.addClass("editablerow");
      var name = td.data("name");
      var type = td.data("type");
      data = { dict_value: name, type: type, value: value };
      console.log(data);
      sendToServer(data);
    }
  });

  const csrftoken = getCookie("csrftoken");
  const sendToServer = (data) => {
    console.log(data);
    $.ajax({
      url: "/dict_update",
      type: "POST",
      headers: { "X-CSRFToken": csrftoken },
      data: data,
    })
      .done((response) => {
        console.log(response);
      })
      .fail(() => {
        console.log("Error occured");
      });
  };
});
