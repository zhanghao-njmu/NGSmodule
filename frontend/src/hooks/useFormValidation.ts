/**
 * Form Validation Hook - Enhanced form validation with real-time feedback
 * 表单验证Hook - 带实时反馈的增强表单验证
 */
import { useState, useCallback } from 'react'
import type { FormInstance } from 'antd'
import { toast } from '../utils/notification'

export interface ValidationRule {
  required?: boolean
  min?: number
  max?: number
  pattern?: RegExp
  validator?: (value: any) => boolean | string
  message?: string
}

export interface ValidationRules {
  [key: string]: ValidationRule | ValidationRule[]
}

export const useFormValidation = (form: FormInstance) => {
  const [errors, setErrors] = useState<Record<string, string>>({})
  const [isDirty, setIsDirty] = useState(false)

  /**
   * Validate a single field
   */
  const validateField = useCallback((fieldName: string, value: any, rules: ValidationRule | ValidationRule[]) => {
    const ruleArray = Array.isArray(rules) ? rules : [rules]

    for (const rule of ruleArray) {
      // Required validation
      if (rule.required && !value) {
        return rule.message || `${fieldName}不能为空`
      }

      // Skip other validations if value is empty and not required
      if (!value) {
        continue
      }

      // Min length validation
      if (rule.min && value.length < rule.min) {
        return rule.message || `${fieldName}至少需要${rule.min}个字符`
      }

      // Max length validation
      if (rule.max && value.length > rule.max) {
        return rule.message || `${fieldName}不能超过${rule.max}个字符`
      }

      // Pattern validation
      if (rule.pattern && !rule.pattern.test(value)) {
        return rule.message || `${fieldName}格式不正确`
      }

      // Custom validator
      if (rule.validator) {
        const result = rule.validator(value)
        if (result !== true) {
          return typeof result === 'string' ? result : rule.message || '验证失败'
        }
      }
    }

    return null
  }, [])

  /**
   * Validate all fields
   */
  const validateForm = useCallback(
    (rules: ValidationRules) => {
      const values = form.getFieldsValue()
      const newErrors: Record<string, string> = {}
      let isValid = true

      Object.keys(rules).forEach((fieldName) => {
        const value = values[fieldName]
        const error = validateField(fieldName, value, rules[fieldName])

        if (error) {
          newErrors[fieldName] = error
          isValid = false
        }
      })

      setErrors(newErrors)

      if (!isValid) {
        toast.warning('请检查表单中的错误')
        // Focus on first error field
        const firstErrorField = Object.keys(newErrors)[0]
        if (firstErrorField) {
          form.scrollToField(firstErrorField)
        }
      }

      return isValid
    },
    [form, validateField],
  )

  /**
   * Handle field change with real-time validation
   */
  const handleFieldChange = useCallback(
    (fieldName: string, value: any, rules: ValidationRule | ValidationRule[]) => {
      setIsDirty(true)

      const error = validateField(fieldName, value, rules)
      setErrors((prev) => {
        const newErrors = { ...prev }
        if (error) {
          newErrors[fieldName] = error
        } else {
          delete newErrors[fieldName]
        }
        return newErrors
      })
    },
    [validateField],
  )

  /**
   * Clear all errors
   */
  const clearErrors = useCallback(() => {
    setErrors({})
    setIsDirty(false)
  }, [])

  /**
   * Check if form has errors
   */
  const hasErrors = Object.keys(errors).length > 0

  return {
    errors,
    isDirty,
    hasErrors,
    validateField,
    validateForm,
    handleFieldChange,
    clearErrors,
  }
}

/**
 * Common validation rules
 */
export const validationRules = {
  required: (message?: string): ValidationRule => ({
    required: true,
    message: message || '此字段为必填项',
  }),

  email: (message?: string): ValidationRule => ({
    pattern: /^[^\s@]+@[^\s@]+\.[^\s@]+$/,
    message: message || '请输入有效的邮箱地址',
  }),

  username: (message?: string): ValidationRule => ({
    pattern: /^[a-zA-Z0-9_-]{3,20}$/,
    message: message || '用户名只能包含字母、数字、下划线和连字符，长度3-20个字符',
  }),

  password: (message?: string): ValidationRule => ({
    min: 8,
    pattern: /^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)[a-zA-Z\d@$!%*?&]{8,}$/,
    message: message || '密码至少8个字符，包含大小写字母和数字',
  }),

  phone: (message?: string): ValidationRule => ({
    pattern: /^1[3-9]\d{9}$/,
    message: message || '请输入有效的手机号码',
  }),

  url: (message?: string): ValidationRule => ({
    pattern: /^https?:\/\/.+/,
    message: message || '请输入有效的URL地址',
  }),

  number: (message?: string): ValidationRule => ({
    pattern: /^\d+$/,
    message: message || '请输入数字',
  }),

  positiveNumber: (message?: string): ValidationRule => ({
    pattern: /^[1-9]\d*$/,
    message: message || '请输入正整数',
  }),

  minLength: (min: number, message?: string): ValidationRule => ({
    min,
    message: message || `至少需要${min}个字符`,
  }),

  maxLength: (max: number, message?: string): ValidationRule => ({
    max,
    message: message || `不能超过${max}个字符`,
  }),

  range: (min: number, max: number, message?: string): ValidationRule => ({
    validator: (value: string) => {
      const num = parseFloat(value)
      return !isNaN(num) && num >= min && num <= max
    },
    message: message || `请输入${min}到${max}之间的数值`,
  }),
}

export default useFormValidation
