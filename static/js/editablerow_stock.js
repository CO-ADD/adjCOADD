$(document).ready(function () {
  $(document).on("dblclick", ".editablerow_stock", function () {
    console.log("editable stock...");
    var value = $(this).text();
    var input =
      "<input type='number' class='input-data' min='0' onkeydown='if(event.keyCode==38) return false;' value='" + value + "' />"+"<button>cancel</button>";
    $(this).html(input);
    
    $(this).removeClass("editablerow_stock");
  });

  $(document).on("blur", ".input-data", function () {
    var value = $(this).val().trim();
    var td = $(this).parent("td");
    $(this).remove();
    td.html(value);
    td.addClass("editablerow_stock");
    var name = td.data("name");
    var type = td.data("type");
    data = { name: name, type: type, value: value };
    console.log(data);
    sendToServer(data, td);
  });
  $(document).on("keypress", ".input-data", function (e) {
    if (e.keyCode === 13) {
      var value = $(this).val().trim();
      var td = $(this).parent("td");
      $(this).remove();
      td.html(value);
      td.addClass("editablerow_stock");
      var name = td.data("name");
      var type = td.data("type");
      data = { name: name, type: type, value: value };
      console.log(data);
      sendToServer(data, td);
    }
  });

  const csrftoken = getCookie("csrftoken");
  const sendToServer = (data, td) => {
    console.log(data);
    $.ajax({
      url: `/organism/updateStock/${data['name']}`,
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
