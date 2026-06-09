# Test helper: return the record with the given `name` from an efa_retention
# object (or NULL if no such record is present).
.retention_record <- function(obj, name) {
  rec <- Filter(function(r) r$name == name, obj$results)
  if (length(rec) == 0) NULL else rec[[1]]
}
